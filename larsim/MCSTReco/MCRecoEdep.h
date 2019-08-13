#ifndef MCRECOEDEP_H
#define MCRECOEDEP_H

// LArSoft
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
namespace sim {
  class SimChannel;
  class SimEnergyDeposit;
}

// Framework includes
namespace fhicl { class ParameterSet; }

// STL
#include <map>
#include <vector>

namespace sim
{

  namespace details {
    // Returns a map with all available plane IDs,
    //  each mapped into an index from a compact range.
    std::map<geo::PlaneID, size_t> createPlaneIndexMap();
  } // namespace details


  class MCEdepHit {

  public:

    //static const unsigned short kINVALID_USHORT;

  public:

    MCEdepHit(){ Clear(); }

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

  class UniquePosition{
  public:
    double _x, _y, _z;

  public:
    UniquePosition(double x=0, double y=0, double z=0)
    { _x = x; _y = y; _z = z; }


    inline bool operator<( const UniquePosition& rhs) const
    {
      if(_x < rhs._x) return true;
      if(rhs._x < _x) return false;
      if(_y < rhs._y) return true;
      if(rhs._y < _y) return false;
      if(_z < rhs._z) return true;
      if(rhs._z < _z) return false;
      return false;
    }

  };


  struct MCEdep {
    struct deposit{
      float energy {};
      float charge {};
      deposit() = default;
      deposit(float e, float c) : energy(e), charge(c) { }
    };

    sim::UniquePosition pos {};
    geo::PlaneID pid {};
    std::vector<deposit> deps {};

    MCEdep() = default;

    MCEdep(sim::UniquePosition p,
           geo::PlaneID pi,
           size_t num_planes,
           float e, float c,
           size_t id) :
           pos(p), pid(pi), deps(num_planes) { deps[id].energy=e; deps[id].charge=c;}
  };

  class MCRecoEdep {

  public:

    /// Default constructor with fhicl parameters
    MCRecoEdep(fhicl::ParameterSet const& pset);
    //ClusterMergeAlg(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    void MakeMCEdep(const std::vector<sim::SimChannel>& schArray);

    void MakeMCEdep(const std::vector<sim::SimEnergyDeposit>& sedArray);

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

    void Clear() {
      _mc_edeps.clear();
      _track_index.clear();
      std::vector<std::vector<sim::MCEdep>>().swap(_mc_edeps);
      std::map<unsigned int,size_t>().swap(_track_index);
  }
  protected:

    std::vector<sim::MCEdep>& __GetEdepArray__(unsigned int track_id);

    bool _debug_mode;
    bool _save_mchit;
    std::map<unsigned int,size_t>      _track_index;
    std::vector<std::vector<sim::MCEdep> > _mc_edeps;

  }; // class MCRecoEdep

} //namespace cluster
#endif
