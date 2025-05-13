#ifndef MCRECOEDEP_H
#define MCRECOEDEP_H

// LArSoft
#include "larcorealg/Geometry/GeometryDataContainers.h" // geo::PlaneDataContainer
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimEnergyDepositLite.h"

// Framework includes
namespace fhicl {
  class ParameterSet;
}

// STL
#include <cassert>
#include <map>
#include <vector>

namespace sim {

  namespace details {
    // Returns a map with all available plane IDs,
    //  each mapped into an index from a compact range.
    geo::PlaneDataContainer<std::size_t> createPlaneIndexMap();
  } // namespace details

  class MCEdepHit {

  public:
    //static const unsigned short kINVALID_USHORT;

  public:
    MCEdepHit() { Clear(); }

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

  struct UniquePosition : geo::Point_t {

    using geo::Point_t::Point_t;
    UniquePosition(geo::Point_t const& p) : geo::Point_t{p} {}

    inline bool operator<(const UniquePosition& rhs) const
    {
      if (X() < rhs.X()) return true;
      if (rhs.X() < X()) return false;
      if (Y() < rhs.Y()) return true;
      if (rhs.Y() < Y()) return false;
      if (Z() < rhs.Z()) return true;
      if (rhs.Z() < Z()) return false;
      return false;
    }
  };

  struct MCEdep {
    struct deposit {
      float energy = 0.0;
      float charge = 0.0;
    };

    sim::UniquePosition pos;
    geo::PlaneID pid;
    std::vector<deposit> deps;

    MCEdep() = default;

    /// @todo Document the protocol for using this object and the (unprotected) `deps` in particular
    MCEdep(sim::UniquePosition p,
           geo::PlaneID const& pi,
           size_t num_planes,
           float e,
           float c,
           size_t id)
      : pos(p), pid(pi), deps(num_planes)
    {
      assert(id < num_planes); // prerequisite
      deps[id].energy = e;
      deps[id].charge = c;
    }
  };

  class MCRecoEdep {

  public:
    /// Default constructor with fhicl parameters
    MCRecoEdep(fhicl::ParameterSet const& pset);

    void MakeMCEdep(const std::vector<sim::SimChannel>& schArray);

    void MakeMCEdep(const std::vector<sim::SimEnergyDeposit>& sedArray);

    void MakeMCEdep(const std::vector<sim::SimEnergyDepositLite>& sedArray);

    bool ExistTrack(const unsigned int track_id) const
    {
      return (_track_index.find(track_id) != _track_index.end());
    }

    /// Converts a track ID to MCEdep array index. Returns -1 if no corresponding array found .
    int TrackToEdepIndex(unsigned int track_id) const
    {
      auto iter = _track_index.find(track_id);
      return (iter == _track_index.end() ? -1 : (int)((*iter).second));
    }

    /// Returns a vector of MCEdep object at the given index
    const std::vector<sim::MCEdep>& GetEdepArrayAt(size_t edep_index) const;

    /// Returns a map of track id <-> MCEdep vector index
    const std::map<unsigned int, size_t> TrackIndexMap() const { return _track_index; }

    void Clear()
    {
      _mc_edeps.clear();
      _track_index.clear();
      std::vector<std::vector<sim::MCEdep>>().swap(_mc_edeps);
      std::map<unsigned int, size_t>().swap(_track_index);
    }

  protected:
    std::vector<sim::MCEdep>& __GetEdepArray__(unsigned int track_id);

    bool _debug_mode;
    bool _save_mchit;
    std::map<unsigned int, size_t> _track_index;
    std::vector<std::vector<sim::MCEdep>> _mc_edeps;

  }; // class MCRecoEdep

} //namespace cluster
#endif
