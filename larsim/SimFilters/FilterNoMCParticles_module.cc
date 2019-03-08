////////////////////////////////////////////////////////////////////////
/// \file  FilterNoMCParticles_module.cc
/// \brief Simple EDFilter filter events with no MCParticles
///
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"

// LArSoft Includes
#include "nusimdata/SimulationBase/MCTruth.h"


namespace simb{
  class MCTruth;
}

///Geant4 interface
namespace simfilter {

  class FilterNoMCParticles : public art::EDFilter
  {
  public:

    explicit FilterNoMCParticles(fhicl::ParameterSet const &pset);

    bool filter(art::Event&) ;

  private:

    std::string fLArG4ModuleLabel;

  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterNoMCParticles::FilterNoMCParticles(fhicl::ParameterSet const& pset) :
    EDFilter{pset},
    fLArG4ModuleLabel    (pset.get< std::string > ("LArG4ModuleLabel"   , "NoLabel")       )
  {}

  //-----------------------------------------------------------------------
  bool FilterNoMCParticles::filter(art::Event& evt)
  {
  //  bool interactionDesired(false);

    art::Handle<std::vector<simb::MCParticle> > mcpHandle;
    evt.getByLabel(fLArG4ModuleLabel,mcpHandle);


    return mcpHandle->size()>0;


  } // end FilterNoMCParticles()function

} // namespace simfilter

namespace simfilter {

  DEFINE_ART_MODULE(FilterNoMCParticles)

} // namespace simfilter
