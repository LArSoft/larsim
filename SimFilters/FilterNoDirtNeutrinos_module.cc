////////////////////////////////////////////////////////////////////////
/// \file  LArG4.cxx
/// \brief Simple EDFilter to require neutrino interaction in TPC
///
/// \version $Id: 
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
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

// LArSoft Includes
#include "MCCheater/BackTracker.h"
#include "Simulation/ParticleList.h"
#include "Simulation/sim.h"
#include "Geometry/Geometry.h"

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
 
  class FilterNoDirtNeutrinos : public art::EDFilter 
  {  
  // explicit EDFilter(ParameterSet const&)  
  public:

    explicit FilterNoDirtNeutrinos(fhicl::ParameterSet const &pset);
    virtual ~FilterNoDirtNeutrinos();                        
    
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
    std::string fGenModuleLabel;
    std::vector<int> fPrimaryVec;
  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterNoDirtNeutrinos::FilterNoDirtNeutrinos(fhicl::ParameterSet const& pset) :
    fG4ModuleLabel    (pset.get< std::string > ("GeantModuleLabel"   , "NoLabel")       )
    , fGenModuleLabel    (pset.get< std::string > ("GenModuleLabel"  , "NoLabel")       )
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Destructor
  FilterNoDirtNeutrinos::~FilterNoDirtNeutrinos() 
  {
  }

  //-----------------------------------------------------------------------
  void FilterNoDirtNeutrinos::beginJob()
  {
    //    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<geo::Geometry> geo;
  
  }

  //-----------------------------------------------------------------------
  void FilterNoDirtNeutrinos::reconfigure(fhicl::ParameterSet const& p)
  {
    fPrimaryVec =  p.get<std::vector<int> >("PrimaryParticles");

    return;
  }

  //-----------------------------------------------------------------------
  bool FilterNoDirtNeutrinos::filter(art::Event& evt) 
  {
    bool interactionDesired(false);
    //get the list of particles from this event
    art::ServiceHandle<cheat::BackTracker> bt;
    art::ServiceHandle<geo::Geometry> geom;

    // * MC truth information
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;

    if (evt.getByLabel(fGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

    for (size_t imc=0; imc<mclist.size() && !interactionDesired ;imc++) 
      {
	bool isMC = !evt.isRealData();
	if (isMC) 
	  { //is MC
	    // GENIE
	    if (!mclist.empty())
	      { //at least one mc record
	
		art::Ptr<simb::MCTruth> mctruth;
		mctruth = mclist.at(imc);

		int nGeniePrimaries(0);
		if (mctruth->NeutrinoSet()) nGeniePrimaries = mctruth->NParticles();
		std::cout << "FilterNoDirtNeutrinos: There are " << nGeniePrimaries << " in the list." ;
	      }
	    // get the particles from the back tracker
	    std::vector<const simb::MCParticle*> pvec;
	    for(size_t i = 0; i < bt->ParticleList().size(); ++i){
	      pvec.push_back(bt->ParticleList().Particle(i));
	    }

	    // Get fiducial volume boundary.
	    double xmin = 0.;
	    double xmax = 2.*geom->DetHalfWidth();
	    double ymin = -geom->DetHalfHeight();
	    double ymax = geom->DetHalfHeight();
	    double zmin = 0.;
	    double zmax = geom->DetLength();
	    

	    for(unsigned int i = 0; i < pvec.size() && !interactionDesired ; ++i)
	      {
		
		const simb::MCParticle* part = pvec.at(i);
		int pdg = abs(part->PdgCode());
		int trackID = part->TrackId();
		
		const art::Ptr<simb::MCTruth> mc = bt->TrackIDToMCTruth(trackID);
		if (mc->Origin() == simb::kBeamNeutrino )
		  {
		    // Now walk through trajectory and see if it enters the TPC
		    int n = part->NumberTrajectoryPoints();
		    for(int j = 0; j < n && !interactionDesired; ++j) 
		      {
			
			TVector3 pos = part->Position(i).Vect();
			if(pos.X() >= xmin &&
			   pos.X() <= xmax &&
			   pos.Y() >= ymin &&
			   pos.Y() <= ymax &&
			   pos.Z() >= zmin &&
			   pos.Z() <= zmax) 
			  {
			    mf::LogInfo("FilterNoDirtNeutrinos") << " Found a Genie-produced particle " << pdg << " in TPC.";
			    interactionDesired = true;
			  }
		      } // end loop on MC trajectory
		  } // end Origin is from beam neutrino, not cosmic.
	      } // end loop over all particles
	  } // end clause that this is true MC
      } // end loop on mclist

    return interactionDesired;
    
  } // end FilterNoDirtNeutrinos()function

} // namespace simfilter

namespace simfilter {

  DEFINE_ART_MODULE(FilterNoDirtNeutrinos)

} // namespace simfilter

#endif // FILTER_FILTERNODIRTNUS_H

