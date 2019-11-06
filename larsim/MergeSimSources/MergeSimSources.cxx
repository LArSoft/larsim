/*!
 * Title:   MergeSimSources Utility Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Class that merges different simulation sources together to created a combined sim list.
 * Typically just merges vectors/maps/etc together. But, if anything as a G4 trackID, applies
 * a user-defined offset to those IDs.
 *
*/

#include <algorithm>
#include <stdexcept>
#include <sstream>

#include "MergeSimSources.h"

sim::MergeSimSourcesUtility::MergeSimSourcesUtility(const std::vector<int>& offsets)
{
  fG4TrackIDOffsets = offsets;
  Reset();
}

void sim::MergeSimSourcesUtility::Reset()
{
  fG4TrackIDRanges.resize(fG4TrackIDOffsets.size(),
			  std::make_pair(std::numeric_limits<int>::max(),
					 std::numeric_limits<int>::min()));
  fMCParticleListMap.resize(fG4TrackIDOffsets.size(),
			    std::vector<size_t>());
}

void sim::MergeSimSourcesUtility::MergeMCParticles( std::vector<simb::MCParticle>& merged_vector,
						    const std::vector<simb::MCParticle>& input_vector,
						    size_t source_index)
{

  if(source_index >= fG4TrackIDOffsets.size())
    std::runtime_error("ERROR in MergeSimSourcesUtility: Source index out of range!");

  fMCParticleListMap[source_index].resize(input_vector.size());
  merged_vector.reserve(merged_vector.size() + input_vector.size());

  std::pair<int,int> range_trackID(std::numeric_limits<int>::max(),
				   std::numeric_limits<int>::min());

  for(size_t i_p=0; i_p<input_vector.size(); i_p++){
    merged_vector.emplace_back(input_vector[i_p],fG4TrackIDOffsets[source_index]);

    fMCParticleListMap[source_index][i_p] = merged_vector.size() - 1;

    if(merged_vector.back().TrackId() < range_trackID.first)
      range_trackID.first = merged_vector.back().TrackId();
    if(merged_vector.back().TrackId() > range_trackID.second)
      range_trackID.second = merged_vector.back().TrackId();

  }

  UpdateG4TrackIDRange(range_trackID,source_index);
}

void sim::MergeSimSourcesUtility::MergeSimChannels(std::vector<sim::SimChannel>& merged_vector,
						   const std::vector<sim::SimChannel>& input_vector,
						   size_t source_index)
{
  if(source_index >= fG4TrackIDOffsets.size())
    std::runtime_error("ERROR in MergeSimSourcesUtility: Source index out of range!");

  merged_vector.reserve( merged_vector.size() + input_vector.size() );

  std::pair<int,int> range_trackID(std::numeric_limits<int>::max(),
				   std::numeric_limits<int>::min());

  for(auto const& simchannel : input_vector){
    std::vector<sim::SimChannel>::iterator it = std::find(merged_vector.begin(),merged_vector.end(),simchannel);

    if(it==merged_vector.end()){
      merged_vector.emplace_back(simchannel.Channel());
      it = merged_vector.end() - 1;
    }

    std::pair<int,int> thisrange = it->MergeSimChannel(simchannel,fG4TrackIDOffsets[source_index]);
    if(thisrange.first < range_trackID.first) range_trackID.first = thisrange.first;
    if(thisrange.second > range_trackID.second) range_trackID.second = thisrange.second;
  }

  UpdateG4TrackIDRange(range_trackID,source_index);
}

void sim::MergeSimSourcesUtility::MergeAuxDetSimChannels(std::vector<sim::AuxDetSimChannel>& merged_vector,
							 const std::vector<sim::AuxDetSimChannel>& input_vector,
							 size_t source_index)
{
  if(source_index >= fG4TrackIDOffsets.size())
    std::runtime_error("ERROR in MergeSimSourcesUtility: Source index out of range!");

  merged_vector.reserve( merged_vector.size() + input_vector.size() );

  std::pair<int,int> range_trackID(std::numeric_limits<int>::max(),
				   std::numeric_limits<int>::min());

  for(auto const& simchannel : input_vector){
    std::vector<sim::AuxDetSimChannel>::iterator it = std::find(merged_vector.begin(),merged_vector.end(),simchannel);

    if(it==merged_vector.end()){
      merged_vector.emplace_back(simchannel.AuxDetID(), simchannel.AuxDetSensitiveID());
      it = merged_vector.end() - 1;
    }

    // re-make the AuxDetSimChannel with both pairs of AuxDetIDEs 
    int offset = fG4TrackIDOffsets[source_index];
    std::vector<sim::AuxDetIDE> all_ides = it->AuxDetIDEs(); 
    for (const sim::AuxDetIDE &ide: simchannel.AuxDetIDEs()) {
      all_ides.emplace_back(ide, offset);

      if( ide.trackID+offset < range_trackID.first  )
        range_trackID.first = ide.trackID+offset;
      if( ide.trackID+offset > range_trackID.second )
        range_trackID.second = ide.trackID+offset;
    }
    

    *it = sim::AuxDetSimChannel(simchannel.AuxDetID(), std::move(all_ides), simchannel.AuxDetSensitiveID());
  }

  UpdateG4TrackIDRange(range_trackID,source_index);
}

void sim::MergeSimSourcesUtility::MergeSimPhotons( std::vector<sim::SimPhotons>& merged_vector,
						   const std::vector<sim::SimPhotons>& input_vector)
{

  merged_vector.reserve( merged_vector.size() + input_vector.size() );

  for(auto const& simphotons : input_vector){
    std::vector<sim::SimPhotons>::iterator it = std::find(merged_vector.begin(),merged_vector.end(),simphotons);

    if(it==merged_vector.end()){
      merged_vector.emplace_back(simphotons.OpChannel());
      it = merged_vector.end() - 1;
    }

    *it += simphotons;
  }
}

void sim::MergeSimSourcesUtility::MergeSimPhotonsLite( std::vector<sim::SimPhotonsLite>& merged_vector,
						       const std::vector<sim::SimPhotonsLite>& input_vector)
{

  merged_vector.reserve( merged_vector.size() + input_vector.size() );

  for(auto const& simphotons : input_vector){
    std::vector<sim::SimPhotonsLite>::iterator it = std::find(merged_vector.begin(),merged_vector.end(),simphotons);

    if(it==merged_vector.end()){
      merged_vector.emplace_back(simphotons.OpChannel);
      it = merged_vector.end() - 1;
    }

    *it += simphotons;
  }
}

void sim::MergeSimSourcesUtility::UpdateG4TrackIDRange(std::pair<int,int> newrange, size_t source_index)
{
  if(source_index >= fG4TrackIDOffsets.size())
    std::runtime_error("ERROR in MergeSimSourcesUtility: Source index out of range!");

  if( newrange.first >= fG4TrackIDRanges[source_index].first &&
      newrange.second <= fG4TrackIDRanges[source_index].second)
    return;

  for(size_t i=0; i<fG4TrackIDRanges.size(); i++){
    if(i==source_index) continue;

    if( (newrange.first >= fG4TrackIDRanges[i].first && newrange.first <= fG4TrackIDRanges[i].second) ||
	(newrange.second >= fG4TrackIDRanges[i].first && newrange.second <= fG4TrackIDRanges[i].second) )
      {
	std::stringstream ss;
	ss << "ERROR in MergeSimSourcesUtility: Source trackIDs overlap!"
	   << "\n\t" << i << "\t" << fG4TrackIDRanges[i].first << " " << fG4TrackIDRanges[i].second
	   << "\n\t" << "n\t" << newrange.first << " " << newrange.second;
	throw std::runtime_error(ss.str());
      }
  }

  if(newrange.first < fG4TrackIDRanges[source_index].first)
    fG4TrackIDRanges[source_index].first = newrange.first;
  if(newrange.second > fG4TrackIDRanges[source_index].second)
    fG4TrackIDRanges[source_index].second = newrange.second;


}
