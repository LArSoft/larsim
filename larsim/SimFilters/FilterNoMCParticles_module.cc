////////////////////////////////////////////////////////////////////////
/// \file  FilterNoMCParticles_module.cc
/// \brief Simple EDFilter filter events with no MCParticles
///
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef FILTER_FILTERNODIRTNUS_H
#define FILTER_FILTERNODIRTNUS_H

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// LArSoft Includes
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/sim.h"
#include "larcore/Geometry/Geometry.h"

// C++ Includes
#include <iostream>
#include <cstring>
#include <sys/stat.h>

namespace simb{
  class MCTruth;
}

namespace sim{
  class ParticleList;
}

///Geant4 interface
namespace simfilter {

  class FilterNoMCParticles : public art::EDFilter
  {
  // explicit EDFilter(ParameterSet const&)
  public:

    explicit FilterNoMCParticles(fhicl::ParameterSet const &pset);
    virtual ~FilterNoMCParticles();

    bool filter(art::Event&) ;
    virtual void reconfigure(fhicl::ParameterSet const&)  ;

    virtual void beginJob()  ;

  private:

    std::string fLArG4ModuleLabel;

  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterNoMCParticles::FilterNoMCParticles(fhicl::ParameterSet const& pset) :
    fLArG4ModuleLabel    (pset.get< std::string > ("LArG4ModuleLabel"   , "NoLabel")       )
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Destructor
  FilterNoMCParticles::~FilterNoMCParticles()
  {
  }

  //-----------------------------------------------------------------------
  void FilterNoMCParticles::beginJob()
  {
    //    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<geo::Geometry> geo;

  }

  //-----------------------------------------------------------------------
  void FilterNoMCParticles::reconfigure(fhicl::ParameterSet const& p)
  {

    return;
  }

  //-----------------------------------------------------------------------
  bool FilterNoMCParticles::filter(art::Event& evt)
  {
    bool interactionDesired(false);

    art::Handle<std::vector<simb::MCParticle> > mcpHandle;
    evt.getByLabel(fLArG4ModuleLabel,mcpHandle);


    return mcpHandle->size()>0;


  } // end FilterNoMCParticles()function

} // namespace simfilter

namespace simfilter {

  DEFINE_ART_MODULE(FilterNoMCParticles)

} // namespace simfilter

#endif // FILTER_FILTERNODIRTNUS_H
