////////////////////////////////////////////////////////////////////////
/// \file  FilterStoppingMuon_module.cc
/// \brief Simple EDFilter to require muon to stop in the TPC
///
/// \author  dcaratelli@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// C++ Includes
#include <iostream>

namespace {

  struct bounds {
    double minBound;
    double maxBound;
  };

  inline bool
  within(bounds const& bnds, double const value)
  {
    return bnds.minBound < value and value < bnds.maxBound;
  }

  class FilterStoppingMuon : public art::SharedFilter {
  public:
    explicit FilterStoppingMuon(fhicl::ParameterSet const& pset,
                                art::ProcessingFrame const&);

  private:
    void beginRun(art::Run&, art::ProcessingFrame const&) override;
    bool filter(art::Event&, art::ProcessingFrame const&) override;
    std::string const fLArG4ModuleLabel;
    bounds fXBounds{};
    bounds fYBounds{};
    bounds fZBounds{};
  };

  //-----------------------------------------------------------------------
  // Constructor
  FilterStoppingMuon::FilterStoppingMuon(fhicl::ParameterSet const& pset,
                                         art::ProcessingFrame const&)
    : SharedFilter{pset}
    , fLArG4ModuleLabel{pset.get<std::string>("LArG4ModuleLabel", "largeant")}
  {
    async<art::InEvent>();
  }

  void
  FilterStoppingMuon::beginRun(art::Run&, art::ProcessingFrame const& frame)
  {
    // Detector geometries are allowed to change on run boundaries.
    auto const geom = frame.serviceHandle<geo::Geometry const>();
    fXBounds = {0., 2. * geom->DetHalfWidth()};
    fYBounds = {-geom->DetHalfHeight(), geom->DetHalfHeight()};
    fZBounds = {0., geom->DetLength()};
  }

  //-----------------------------------------------------------------------
  bool
  FilterStoppingMuon::filter(art::Event& evt, art::ProcessingFrame const&)
  {
    // get the particles produced by largeant
    auto const& mcps =
      *evt.getValidHandle<std::vector<simb::MCParticle>>(fLArG4ModuleLabel);

    for (auto const& part : mcps) {
      // skip anything that isn't a muon
      if (std::abs(part.PdgCode()) != 13) {
        continue;
      }

      // make sure the end position is inside the TPC
      if (within(fXBounds, part.EndX()) and within(fYBounds, part.EndY()) and
          within(fZBounds, part.EndZ())) {
        // we found a stopping muon -> return true (keep this event)
        // N.B. Printing this out during multi-threaded execution may
        // result in mangled printouts.  If that is a problem, then
        // these printouts should be made only if a 'debug' mode is
        // specified, in which case 'serialize()' would be called
        // instead of 'async'.
        std::cout << "************* IN TPC *******************" << std::endl;
        return true;
      }

    } // for all mcparticles

    return false;

  } // end FilterStoppingMuon()function

}

DEFINE_ART_MODULE(FilterStoppingMuon)
