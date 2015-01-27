////////////////////////////////////////////////////////////////////////
//
//  MCShowerRecoPart source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCSHOWERRECOPART_CXX
#define MCSHOWERRECOPART_CXX

#include "MCShowerRecoPart.h"

namespace sim {

  const unsigned int MCShowerRecoPart::kINVALID_UINT = std::numeric_limits<unsigned int>::max();
  const int MCShowerRecoPart::kINVALID_INT = std::numeric_limits<int>::max();

  //##################################################################
  MCShowerRecoPart::MCShowerRecoPart(fhicl::ParameterSet const& pset)
  //##################################################################
  {
    _debug_mode = pset.get<bool>("DebugMode");
  }

  //----------------------------------------------------------------------------------------------- 
  void MCShowerRecoPart::ConstructShower(const MCRecoPart& part_v)
  //----------------------------------------------------------------------------------------------- 
  {
    if(!part_v.size()) return;

    _shower_id.clear();
    _shower_id.resize(part_v.size(),-1);
    _shower_index.clear();
    _shower_daughters.clear();

    // Construct MCShower
    std::vector<std::multimap<double,unsigned int> > daughter_map;
    for(size_t i=0; i<part_v.size(); ++i) {

      auto const& mcp = part_v[i];

      int candidate_mom_index=-1;
      if( mcp._pdgcode == 22 ||
          mcp._pdgcode == 11 ||
          mcp._pdgcode == -11 )
	candidate_mom_index = i;

      unsigned int mom_track = mcp._mother;
      auto mom_iter = part_v._track_index.find(mom_track);
      while(mom_iter != part_v._track_index.end()) {

        unsigned int mom_index = (*mom_iter).second;

        if( part_v.at(mom_index)._pdgcode == 22 || part_v.at(mom_index)._pdgcode == 11 || part_v.at(mom_index)._pdgcode == -11 )

          candidate_mom_index = mom_index;
	
        mom_iter = part_v._track_index.find(part_v.at(mom_index)._mother);

      }

      if(candidate_mom_index >= 0) {
	
	auto candidate_mom_iter = _shower_index.find(candidate_mom_index);
	if(candidate_mom_iter == _shower_index.end()) {
	  _shower_index.insert(std::make_pair((unsigned int)candidate_mom_index, (unsigned int)_shower_index.size()));
	  daughter_map.push_back(std::multimap<double,unsigned int>());
	}
	unsigned int shower_index = (*_shower_index.find(candidate_mom_index)).second;
	daughter_map.at(shower_index).insert(std::make_pair((double)(mcp._start_vtx[3]),(unsigned int)i));
	_shower_id.at(i) = shower_index;
	
      } else if(_debug_mode) {

	std::cout
	  << "Found a particle that does not belong to a shower!" << std::endl
	  << Form(" PDGID: %d ... Track %d @ (%g,%g,%g,%g) with (%g,%g,%g,%g)",
		  mcp._pdgcode,
		  mcp._track_id,
		  mcp._start_vtx[0],
		  mcp._start_vtx[1],
		  mcp._start_vtx[2],
		  mcp._start_vtx[3],
		  mcp._start_mom[0],
		  mcp._start_mom[1],
		  mcp._start_mom[2],
		  mcp._start_mom[3])
	  << std::endl << std::endl;

      }
    }


    if(_debug_mode)
      std::cout
	<< Form("Found %zu MCShowers....",_shower_index.size()) << std::endl;

    _shower_daughters.resize(_shower_index.size(),std::vector<unsigned int>());
    for(const auto &mom : _shower_index) {

      _shower_daughters.at(mom.second).reserve(daughter_map.at(mom.second).size());
      for(auto const &part_index : daughter_map.at(mom.second))

	_shower_daughters.at(mom.second).push_back(part_index.second);

      auto const& mcp = part_v.at(mom.first);
      if(_debug_mode) 
	std::cout 
	  << Form("PDGID: %d ... Track %d @ (%g,%g,%g,%g) with (%g,%g,%g,%g) ... %zu daughters!",
		  mcp._pdgcode,
		  mcp._track_id,
		  mcp._start_vtx[0],
		  mcp._start_vtx[1],
		  mcp._start_vtx[2],
		  mcp._start_vtx[3],
		  mcp._start_mom[0],
		  mcp._start_mom[1],
		  mcp._start_mom[2],
		  mcp._start_mom[3],
		  _shower_daughters.at(mom.second).size())
	  << std::endl;
    }

    if(_debug_mode)
      std::cout<<std::endl;
  }
  
}
#endif
