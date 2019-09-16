////////////////////////////////////////////////////////////////////////
/// \file  FilterNoMCParticles_module.cc
/// \brief Simple EDFilter filter events with no MCParticles
///
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft Includes
#include "nusimdata/SimulationBase/MCParticle.h"

#include <vector>

/// Geant4 interface
namespace simfilter {

  class FilterNoMCParticles : public art::SharedFilter {
  public:
    explicit FilterNoMCParticles(fhicl::ParameterSet const& pset,
                                 art::ProcessingFrame const&);

  private:
    bool filter(art::Event&, art::ProcessingFrame const&) override;
    std::string const fLArG4ModuleLabel;
  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterNoMCParticles::FilterNoMCParticles(fhicl::ParameterSet const& pset,
                                           art::ProcessingFrame const&)
    : SharedFilter{pset}
    , fLArG4ModuleLabel{pset.get<std::string>("LArG4ModuleLabel", "NoLabel")}
  {
    async<art::InEvent>();
  }

  //-----------------------------------------------------------------------
  bool
  FilterNoMCParticles::filter(art::Event& evt, art::ProcessingFrame const&)
  {
    auto const& mcps =
      *evt.getValidHandle<std::vector<simb::MCParticle>>(fLArG4ModuleLabel);
    return not mcps.empty();
  }

} // namespace simfilter

DEFINE_ART_MODULE(simfilter::FilterNoMCParticles)
