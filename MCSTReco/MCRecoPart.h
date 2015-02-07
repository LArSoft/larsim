#ifndef MCRECOPART_H
#define MCRECOPART_H

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// LArSoft
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Geometry/Geometry.h"
#include "MCBase/MCLimits.h"

// STL
#include <set>
#include <vector>
#include <sstream>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace sim
{

  class MCMiniPart {

  public:

    MCMiniPart() {Reset();}

    virtual ~MCMiniPart(){}

    unsigned int   _track_id;
    std::string    _process;
    unsigned int   _mother;
    unsigned int   _ancestor;
    int            _pdgcode;
    TLorentzVector _start_vtx;
    TLorentzVector _start_mom;
    TLorentzVector _end_vtx;
    TLorentzVector _end_mom;
    std::vector<std::pair<TLorentzVector,TLorentzVector> > _det_path;
    std::set<unsigned int> _daughters;
    ::simb::Origin_t _origin;

    void Reset(){
      _track_id = _mother = _ancestor = kINVALID_UINT;
      _pdgcode  = kINVALID_INT;
      _process  = "";
      _origin   = ::simb::kUnknown;

      TLorentzVector invalid(kINVALID_DOUBLE,
			     kINVALID_DOUBLE,
			     kINVALID_DOUBLE,
			     kINVALID_DOUBLE);
      _start_vtx = invalid;
      _start_mom = invalid;
      _end_vtx = invalid;
      _end_mom = invalid;
      _daughters.clear();
      _det_path.clear();
    }

  };

  class MCRecoPart : public std::vector<sim::MCMiniPart> {

  public:

    /// Default constructor with fhicl parameters
    MCRecoPart(fhicl::ParameterSet const& pset);

    /// Default destructor
    virtual ~MCRecoPart(){};

    void AddParticles(const std::vector<simb::MCParticle>& mcp_v,
		      const std::vector<simb::Origin_t>&   orig_v);

    unsigned int AncestorTrackID(const unsigned int part_index) const;

    unsigned int MotherTrackID(const unsigned int part_index) const;

    /*
      Take TrackID and returns the corresponding particle unique index number (MCParticle array index)
      Returns kINVALID_UINT if nothing found. 
    */
    unsigned int TrackToParticleIndex(const unsigned int track_id) const
    { 
      auto const iter (_track_index.find(track_id));
      if(iter==_track_index.end()) return kINVALID_UINT;
      return (*iter).second;
    }

    bool InDetector(const double& x,
		    const double& y,
		    const double& z) const;
    
  public:

    /// Track ID => Index Map
    std::map<unsigned int, unsigned int> _track_index;
    //std::map<size_t,size_t> _track_index;

    /// PDG code list for which particle's trajectory within the detector is saved
    std::set<int> _pdg_list;

  protected:

    double _x_max, _x_min;
    double _y_max, _y_min;
    double _z_max, _z_min;

  }; // class MCRecoPart
  
} //namespace cluster
#endif
