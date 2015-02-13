#ifndef MERGESIMSOURCES_H
#define MERGESIMSOURCES_H

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

#include "Simulation/SimChannel.h"
#include "Simulation/SimPhotons.h"
#include "Simulation/AuxDetSimChannel.h"

namespace sim{

  class MergeSimSourcesUtility{

  public:
    
    void MergeSimChannels( std::vector<sim::SimChannel>&,
			   const std::vector<sim::SimChannel>&,
			   size_t);

    void MergeAuxDetSimChannels( std::vector<sim::AuxDetSimChannel>&,
				 const std::vector<sim::AuxDetSimChannel>&,
				 size_t);

    void MergeSimPhotons( std::vector<sim::SimPhotons>&,
			  const std::vector<sim::SimPhotons>&);

    void MergeSimPhotonsLite( std::vector<sim::SimPhotonsLite>&,
			      const std::vector<sim::SimPhotonsLite>&);

  private:

    std::vector<int>                   fG4TrackIDOffsets;
    std::vector< std::pair<int,int> >  fG4TrackIDRanges;

    void UpdateG4TrackIDRange(std::pair<int,int>,size_t);
    
  }; //end MergeSimSourcesUtility class

} //end namespace sim

#endif
