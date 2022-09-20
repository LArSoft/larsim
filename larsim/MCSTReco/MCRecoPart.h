#ifndef MCRECOPART_H
#define MCRECOPART_H

// ART includes
namespace fhicl { class ParameterSet; }

// LArSoft
#include "lardataobj/MCBase/MCLimits.h" // kINVALID_X
#include "lardataobj/MCBase/MCMiniPart.h" // sim::MCMiniPart
#include "lardataobj/MCBase/MCParticleLite.h" // sim::MCParticleLite
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h" // simb::Origin_t

// STL
#include <set>
#include <utility> // std::pair<>
#include <vector>

#include "TLorentzVector.h"

namespace sim
{
  class MCRecoPart : public std::vector<sim::MCMiniPart> {

  public:

    /// Default constructor with fhicl parameters
    MCRecoPart(fhicl::ParameterSet const& pset);

    /// Default destructor
    virtual ~MCRecoPart(){};

    void AddParticles(const std::vector<simb::MCParticle>& mcp_v,
                      const std::vector<simb::Origin_t>&   orig_v,
                      const std::vector<sim::MCParticleLite>&  mcmp_v = {});

    unsigned int AncestorTrackID(const unsigned int part_index);

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

    double _x_max; //!< x-max of volume box used to determine whether to save track information
    double _x_min; //!< x-min of volume box used to determine whether to save track information
    double _y_max; //!< y-max of volume box used to determine whether to save track information
    double _y_min; //!< y-min of volume box used to determine whether to save track information
    double _z_max; //!< z-max of volume box used to determine whether to save track information
    double _z_min; //!< z-min of volume box used to determine whether to save track information

  }; // class MCRecoPart

} //namespace cluster
#endif
