////////////////////////////////////////////////////////////////////////
/// \file  LArG4.cxx
/// \brief Simple EDFilter to require a particular pdg is present as a primary
///
/// \version $Id: 
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef FILTER_FILTERPRIMARYPDG_H
#define FILTER_FILTERPRIMARYPDG_H 

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
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

// LArSoft Includes
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/Simulation/ParticleList.h"
#include "larsim/Simulation/sim.h"
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
 
  class FilterPrimaryPDG : public art::EDFilter 
  {  
  // explicit EDFilter(ParameterSet const&)  
  public:

    explicit FilterPrimaryPDG(fhicl::ParameterSet const &pset);
    virtual ~FilterPrimaryPDG();                        
    
    bool filter(art::Event&) ;
    virtual void reconfigure(fhicl::ParameterSet const&)  ;
      
    virtual void beginJob()  ;
    /*
    virtual void endJob()  ;
    virtual bool beginRun(art::Run &)  ;
    virtual bool endRun(art::Run &)  ;
    virtual bool beginSubRun(art::SubRun &)  ;
    virtual bool endSubRun(art::SubRun &)  ;
    */
    private:

    std::string fG4ModuleLabel;
    std::vector<int> fPrimaryVec;
    /*
  virtual void respondToOpenInputFile(FileBlock const& fb)  
  virtual void respondToCloseInputFile(FileBlock const& fb)  
  virtual void respondToOpenOutputFiles(FileBlock const& fb)  
  virtual void respondToCloseOutputFiles(FileBlock const& fb)  
    */
  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterPrimaryPDG::FilterPrimaryPDG(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Destructor
  FilterPrimaryPDG::~FilterPrimaryPDG() 
  {
  }

  //-----------------------------------------------------------------------
  void FilterPrimaryPDG::beginJob()
  {
    //    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<geo::Geometry> geo;
  
  }

  //-----------------------------------------------------------------------
  void FilterPrimaryPDG::reconfigure(fhicl::ParameterSet const& p)
  {
    fPrimaryVec =  p.get<std::vector<int> >("PrimaryParticles");
    return;
  }

  //-----------------------------------------------------------------------
  bool FilterPrimaryPDG::filter(art::Event& evt) 
  {

    //get the list of particles from this event
    art::ServiceHandle<cheat::BackTracker> bt;
    art::ServiceHandle<geo::Geometry> geom;

    // get the particles from the back tracker
    const sim::ParticleList& Particles = bt->ParticleList();
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

#endif // FILTER_FILTERPRIMARYPDG_H

