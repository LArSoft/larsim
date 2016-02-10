////////////////////////////////////////////////////////////////////////
/// \file  LArG4.cxx
/// \brief Simple EDFilter to require muon to stop in the TPC
///
/// \version $Id: 
/// \author  dcaratelli@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
#ifndef FILTER_FILTERSTOPPINGMUON_H
#define FILTER_FILTERSTOPPINGMUON_H 

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
 
  class FilterStoppingMuon : public art::EDFilter 
  {  
  // explicit EDFilter(ParameterSet const&)  
  public:

    explicit FilterStoppingMuon(fhicl::ParameterSet const &pset);
    virtual ~FilterStoppingMuon();                        
    
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
  FilterStoppingMuon::FilterStoppingMuon(fhicl::ParameterSet const& pset) :
    fLArG4ModuleLabel    (pset.get< std::string > ("LArG4ModuleLabel"   , "largeant")       )
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Destructor
  FilterStoppingMuon::~FilterStoppingMuon() 
  {
  }

  //-----------------------------------------------------------------------
  void FilterStoppingMuon::beginJob()
  {
    //    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<geo::Geometry> geo;
  
  }

  //-----------------------------------------------------------------------
  void FilterStoppingMuon::reconfigure(fhicl::ParameterSet const& p)
  {
    
    return;
  }

  //-----------------------------------------------------------------------
  bool FilterStoppingMuon::filter(art::Event& evt) 
  {
    bool interactionDesired(false);
    //get the list of particles from this event
    art::ServiceHandle<geo::Geometry> geom;

    // * MC truth information
    art::Handle<std::vector<simb::MCParticle> > mcpHandle;

    // get the particles produced by largeant
    evt.getByLabel(fLArG4ModuleLabel,mcpHandle);
    
    double xmin = 0.;
    double xmax = 2.*geom->DetHalfWidth();
    double ymin = -geom->DetHalfHeight();
    double ymax = geom->DetHalfHeight();
    double zmin = 0.;
    double zmax = geom->DetLength();
    
    for(size_t i=0; i < mcpHandle->size(); ++i) {
      
      const simb::MCParticle* part(&mcpHandle->at(i));
      int pdg = part->PdgCode();
      
      // skip anything that isn't a muon
      if ( (pdg != 13) and (pdg != -13) )
	continue;
      
      // get the end position
      double endX = part->EndX();
      double endY = part->EndY();
      double endZ = part->EndZ();
      
      // make sure the end is inside the TPC
      if ( (endX > xmin) and (endX < xmax) and
	   (endY > ymin) and (endY < ymax) and
	   (endZ > zmin) and (endZ < zmax) ){
	// we found a stopping muon -> return true (keep this event)
	std::cout << "************* IN TPC *******************" << std::endl;
	return true;
      }
      
    }// for all mcparticles
    
    return false;
    
  } // end FilterStoppingMuon()function
  
} // namespace simfilter

namespace simfilter {
  
  DEFINE_ART_MODULE(FilterStoppingMuon)

} // namespace simfilter

#endif // FILTER_FILTERNODIRTNUS_H

