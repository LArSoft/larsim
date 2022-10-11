#ifndef MCRECOPART_H
#define MCRECOPART_H

// ART includes
namespace fhicl {
  class ParameterSet;
}

// LArSoft
#include "lardataobj/MCBase/MCLimits.h"       // kINVALID_X
#include "lardataobj/MCBase/MCParticleLite.h" // sim::MCParticleLite
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h" // simb::Origin_t

// STL
#include <set>
#include <utility> // std::pair<>
#include <vector>

#include "TLorentzVector.h"

namespace sim {

  class MCMiniPart {

  public:
    MCMiniPart() { Reset(); }

    virtual ~MCMiniPart() {}

    unsigned int _track_id;
    std::string _process;
    unsigned int _mother;
    unsigned int _ancestor;
    int _pdgcode;
    TLorentzVector _start_vtx;
    TLorentzVector _start_mom;
    TLorentzVector _end_vtx;
    TLorentzVector _end_mom;
    std::vector<std::pair<TLorentzVector, TLorentzVector>> _det_path;
    std::set<unsigned int> _daughters;
    ::simb::Origin_t _origin;

    void Reset()
    {
      _track_id = _mother = _ancestor = kINVALID_UINT;
      _pdgcode = kINVALID_INT;
      _process = "";
      _origin = ::simb::kUnknown;

      TLorentzVector invalid(kINVALID_DOUBLE, kINVALID_DOUBLE, kINVALID_DOUBLE, kINVALID_DOUBLE);
      _start_vtx = invalid;
      _start_mom = invalid;
      _end_vtx = invalid;
      _end_mom = invalid;
      _daughters.clear();
      _det_path.clear();
    }

    bool HasDaughter(unsigned int d) const
    {
      return std::find(_daughters.begin(), _daughters.end(), d) != _daughters.end();
    }

    void AddDaughter(unsigned int d) { _daughters.insert(d); }

    operator simb::MCParticle() const
    {
      simb::MCParticle mcpart(_track_id, _pdgcode, _process, _mother);
      // Also converting back from LArSoft units
      mcpart.AddTrajectoryPoint(_start_vtx, 1.e-3 * _start_mom);
      mcpart.AddTrajectoryPoint(_end_vtx, 1.e-3 * _end_mom);
      for (auto const& d : _daughters)
        mcpart.AddDaughter(d);
      return mcpart;
    }

    MCMiniPart(const simb::MCParticle& p)
    {
      Reset();
      _track_id = p.TrackId();
      _pdgcode = p.PdgCode();
      _mother = p.Mother();
      _process = p.Process();
      _start_vtx = p.Position();
      _start_mom = 1.e3 * p.Momentum(); // Change units to (MeV, cm, us)
      _end_vtx = p.EndPosition();
      _end_mom = 1.e3 * p.EndMomentum(); // idem as above
    }

    MCMiniPart(const sim::MCParticleLite& p)
    {
      Reset();
      _track_id = p.TrackID();
      _pdgcode = p.PdgCode();
      _mother = p.Mother();
      _process = p.Process();
      _start_vtx = p.StartVtx();
      _start_mom = 1.e3 * p.StartMom(); // Change units to (MeV, cm, us)
      _end_vtx = p.EndVtx();
      _end_mom = 1.e3 * p.EndMom(); // idem as above
    }
  };

  class MCRecoPart : public std::vector<sim::MCMiniPart> {

  public:
    /// Default constructor with fhicl parameters
    MCRecoPart(fhicl::ParameterSet const& pset);

    /// Default destructor
    virtual ~MCRecoPart(){};

    void AddParticles(const std::vector<simb::MCParticle>& mcp_v,
                      const std::vector<simb::Origin_t>& orig_v,
                      const std::vector<sim::MCParticleLite>& mcmp_v = {});

    unsigned int AncestorTrackID(const unsigned int part_index);

    unsigned int MotherTrackID(const unsigned int part_index) const;

    /*
      Take TrackID and returns the corresponding particle unique index number (MCParticle array index)
      Returns kINVALID_UINT if nothing found.
    */
    unsigned int TrackToParticleIndex(const unsigned int track_id) const
    {
      auto const iter(_track_index.find(track_id));
      if (iter == _track_index.end()) return kINVALID_UINT;
      return (*iter).second;
    }

    bool InDetector(const double& x, const double& y, const double& z) const;

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
