////////////////////////////////////////////////////////////////////////
/// \file  FilterPrimaryPDG_module.cc
/// \brief Simple EDFilter to require a particular pdg is present as a primary
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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// LArSoft Includes
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "lardataobj/Simulation/sim.h"
#include "larcore/Geometry/Geometry.h"

// C++ Includes
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

  class FilterPrimaryPDG : public art::EDFilter
  {
  public:

    explicit FilterPrimaryPDG(fhicl::ParameterSet const &pset);

    private:
    bool filter(art::Event&) override;

    std::string fG4ModuleLabel;
    std::vector<int> fPrimaryVec;
  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterPrimaryPDG::FilterPrimaryPDG(fhicl::ParameterSet const& pset)
    : EDFilter{pset}, fPrimaryVec{pset.get<std::vector<int> >("PrimaryParticles")}
  {}

  //-----------------------------------------------------------------------
  bool FilterPrimaryPDG::filter(art::Event& evt)
  {

    //get the list of particles from this event
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    art::ServiceHandle<geo::Geometry const> geom;

    // get the particles from the back tracker
    const sim::ParticleList& Particles = pi_serv->ParticleList();
    std::vector<const simb::MCParticle*> pvec;
    pvec.reserve(Particles.size());
    for (const auto& PartPair: Particles) {
      pvec.push_back(PartPair.second);
      // fPDGCodes->Fill(PartPair.second->PdgCode());
    }

    bool pdgDesired(false);
    for(unsigned int i = 0; i < pvec.size(); ++i)
      {
	for (int pdg : fPrimaryVec)
	  {
	    const std::string sprim("primary");
	    if(pvec[i]->PdgCode() == pdg)
	      {
		Char_t tProcess[50];
		for(unsigned int s = 0; s < pvec[i]->Process().length(); ++s)
		  *(tProcess+s) = pvec[i]->Process()[s];
		std::string sProcess(tProcess);
		if (!sProcess.compare(sprim))
		  {
		    mf::LogInfo("FilterPrimaryPDG") << " Found a primary " << pdg << " in event.";
		    pdgDesired = true;
		  }
	      }
	  }

      }

    return pdgDesired;
  }

} // namespace simfilter

namespace simfilter {

  DEFINE_ART_MODULE(FilterPrimaryPDG)

} // namespace simfilter
