////////////////////////////////////////////////////////////////////////
/// \file  FilterCryostatNus_module.cc
/// \brief EDFilter to require projected generator trajectories in volumes
/// within a particular time window.
///
/// \author  Matthew.Bass@physics.ox.ac.uk
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// C++ Includes
#include <cmath> // std::abs()

namespace simfilter {

  class FilterCryostatNus : public art::EDFilter {
  public:
    explicit FilterCryostatNus(fhicl::ParameterSet const& pset);

  private:
    bool filter(art::Event&) override;
    bool fKeepNusInCryostat; // true: keep cryostat nuint; false: filter them
  };

} // namespace simfilter

namespace simfilter {

  FilterCryostatNus::FilterCryostatNus(fhicl::ParameterSet const& pset)
    : EDFilter{pset}
    , fKeepNusInCryostat{pset.get<bool>("KeepNusInCryostat", false)}
  {}

  bool
  FilterCryostatNus::filter(art::Event& evt)
  {
    // get the list of particles from this event
    auto const& geom = *(lar::providerFrom<geo::Geometry>());

    std::vector<art::Handle<std::vector<simb::MCTruth>>> allmclists;
    evt.getManyByType(allmclists);

    bool inCryostatNu = false;
    for (auto const& mclistHandle : allmclists) {
      for (simb::MCTruth const& mct : *mclistHandle) {

        // get nu, does it end in cyrostat?
        for (int ipart = 0; ipart < mct.NParticles(); ipart++) {
          auto const& part = mct.GetParticle(ipart);
          auto const absPDGID = std::abs(part.PdgCode());
          if (absPDGID == 12 || absPDGID == 14 || absPDGID == 16) {
            const TLorentzVector& end4 = part.EndPosition();
            if (geom.PositionToCryostatPtr({end4.X(), end4.Y(), end4.Z()}) !=
                nullptr) {
              inCryostatNu = true;
            }
          } // if neutrino
        }   // for particles

      } // end loop over mctruth col

    } // end loop over all mctruth lists

    return fKeepNusInCryostat ^ (!inCryostatNu);
  }

} // namespace simfilter

DEFINE_ART_MODULE(simfilter::FilterCryostatNus)
