#ifndef MCRECOEDEP_H
#define MCRECOEDEP_H

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
#include "Simulation/SimChannel.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Geometry/Geometry.h"
// STL
#include <set>
#include <vector>
#include <sstream>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace sim
{

  class MCEdepHit {
    
  public:

    //static const unsigned short kINVALID_USHORT;

  public:

    MCEdepHit(){ Clear(); }
    ~MCEdepHit(){}
    
    unsigned short timeStart;
    unsigned short timeEnd;
    unsigned short timeMax;
    float qSum;
    float qMax;
    void Clear()
    {
      //timeStart = timeEnd = timeMax = kINVALID_USHORT;
      timeStart = timeEnd = timeMax = 0;
      qSum = qMax = 0;
    }
  };

  class MCEdep {

  public:
    //static const short kINVALID_SHORT;
    
  public:
    MCEdep(){ Clear(); }
    ~MCEdep(){}
    /*
    short x;
    short y;
    short z;
    */
    float x;
    float y;
    float z;
    float energy;
    double qU;
    double qV;
    double qW;
    std::map<unsigned short,sim::MCEdepHit> mchits;
    void Clear() {
      //x=y=z=kINVALID_SHORT;
      x=y=z=0;
      energy=0;
      qU=qV=qW=0;
      mchits.clear();
    }

  };

  class MCRecoEdep {

  public:

    /// Default constructor with fhicl parameters
    MCRecoEdep(fhicl::ParameterSet const& pset);
    //ClusterMergeAlg(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    /// Default destructor
    virtual ~MCRecoEdep(){};

    void MakeMCEdep(const std::vector<sim::SimChannel>& schArray);

    bool ExistTrack(const unsigned int track_id) const 
    { return (_track_index.find(track_id) != _track_index.end()); } 

    /// Converts a track ID to MCEdep array index. Returns -1 if no corresponding array found .
    int TrackToEdepIndex(unsigned int track_id) const
    { 
      auto iter = _track_index.find(track_id);
      return (iter == _track_index.end() ? -1 : (int)((*iter).second));
    }

    /// Returns a vector of MCEdep object at the given index
    const std::vector<sim::MCEdep>& GetEdepArrayAt(size_t edep_index) const;

    /// Returns a map of track id <-> MCEdep vector index
    const std::map<unsigned int,size_t> TrackIndexMap() const
    { return _track_index; }

  protected:

    std::vector<sim::MCEdep>& __GetEdepArray__(unsigned int track_id);

    bool _debug_mode;
    bool _save_mchit;
    std::map<unsigned int,size_t>      _track_index;
    std::vector<std::vector<sim::MCEdep> > _mc_edeps;
    
  }; // class MCRecoEdep
  
} //namespace cluster
#endif
