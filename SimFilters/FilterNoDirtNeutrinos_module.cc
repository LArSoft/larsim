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
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/FindOneP.h"

// LArSoft Includes
#include "Simulation/ParticleList.h"
#include "SimulationBase/MCTruth.h"
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

    std::string fLArG4ModuleLabel;
    std::string fGenModuleLabel;
    bool        fKeepCryostatNeutrinos;

  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterNoDirtNeutrinos::FilterNoDirtNeutrinos(fhicl::ParameterSet const& pset) :
    fLArG4ModuleLabel    (pset.get< std::string > ("LArG4ModuleLabel"   , "NoLabel")       )
    , fGenModuleLabel    (pset.get< std::string > ("GenModuleLabel"  , "NoLabel")       )
    , fKeepCryostatNeutrinos    (pset.get< bool > ("KeepCryostatNeutrinos", false)      )
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

    return;
  }

  //-----------------------------------------------------------------------
  bool FilterNoDirtNeutrinos::filter(art::Event& evt) 
  {
    bool interactionDesired(false);
    //get the list of particles from this event
    art::ServiceHandle<geo::Geometry> geom;

    // * MC truth information
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    art::Handle<std::vector<simb::MCParticle> > mcpHandle;



    if (evt.getByLabel(fGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);
    evt.getByLabel(fLArG4ModuleLabel,mcpHandle);
    
    //    std::cout << "FilterNoDirtNeutrinos: mclist.size() is " << mclist.size()<< std::endl ;
    
  std::set<art::Ptr<simb::MCTruth> > mctSetGENIE;
  for(size_t i=0; i<mctruthListHandle->size(); ++i) 
    {
      art::Ptr<simb::MCTruth> mct_ptr(mctruthListHandle,i);
      if( mctSetGENIE.find(mct_ptr) == mctSetGENIE.end() ) mctSetGENIE.insert(mct_ptr);
    }

  // Get the MCTruths from associations to our particles
  art::FindOneP<simb::MCTruth> assMCT(mcpHandle, evt, "largeant");

  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;

  if(fKeepCryostatNeutrinos)
  {
    // Get cryostat (box) volume boundary.
    xmin = geom->DetHalfWidth() - geom->CryostatHalfWidth();
    xmax = geom->DetHalfWidth() + geom->CryostatHalfWidth();
    ymin = -geom->CryostatHalfHeight();
    ymax = geom->CryostatHalfHeight();
    zmin = geom->DetLength()/2. - geom->DetLength()/2.;
    zmax = geom->DetLength()/2. + geom->DetLength()/2.;
  }
  else
  {
    // Get fiducial volume boundary.
    xmin = 0.;
    xmax = 2.*geom->DetHalfWidth();
    ymin = -geom->DetHalfHeight();
    ymax = geom->DetHalfHeight();
    zmin = 0.;
    zmax = geom->DetLength();
  }

  //  std::cout << "FilterNoDirtNeutrinos: mcpHandle->size() is " << mcpHandle->size()<< std::endl ;
  // Now let's loop over G4 MCParticle list and track back MCTruth    
  bool inTPC (false);
  for(size_t i=0; i < mcpHandle->size() && !inTPC; ++i) 
    {
      const art::Ptr<simb::MCParticle> mcp_ptr(mcpHandle,i);
      const art::Ptr<simb::MCTruth> &mct = assMCT.at(i);
      if( mctSetGENIE.find(mct) == mctSetGENIE.end() ) 
	{
	  // This is non-genie
	  continue;
	}
      else
	{
	// This is genie

	  const simb::MCParticle* part(&mcpHandle->at(i));
	  int pdg = part->PdgCode();
	  int trackID = part->TrackId();

	  //	std::cout << "FilterNoDirtNeutrinos: i is " << i << std::endl ;
	  // Now walk through trajectory and see if it enters the TPC
	  int n = part->NumberTrajectoryPoints();
	  for(int j = 0; j < n && !inTPC; ++j) 
	    {
	      //	    std::cout << "FilterNoDirtNeutrinos: Loop  counter on NumTrajPt j is " << j << std::endl ;		
	      
	      TVector3 pos = part->Position(j).Vect();
	      if(pos.X() >= xmin &&
		 pos.X() <= xmax &&
		 pos.Y() >= ymin &&
		 pos.Y() <= ymax &&
		 pos.Z() >= zmin &&
		 pos.Z() <= zmax) 
		{
		  interactionDesired = true;
		  //		  std::cout << "FilterNoDirtNeutrinos: Genie daughter found in TPC. G4Particle " << i << " , TrackID/pdg " << trackID << "/ " << pdg << " is discovered." << std::endl ;		
		  std::cout << "FilterNoDirtNeutrinos: Genie daughter found in TPC. G4Particle " << std::endl ;		
		  inTPC=true;
		}
	    } // trajectory loop
	} // end Genie particle
    } // loop on MCPHandle

  return interactionDesired;
    
  } // end FilterNoDirtNeutrinos()function

} // namespace simfilter

namespace simfilter {

  DEFINE_ART_MODULE(FilterNoDirtNeutrinos)

} // namespace simfilter

#endif // FILTER_FILTERNODIRTNUS_H

