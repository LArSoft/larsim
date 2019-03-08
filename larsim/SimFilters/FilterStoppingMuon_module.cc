////////////////////////////////////////////////////////////////////////
/// \file  FilterStoppingMuon_module.cc
/// \brief Simple EDFilter to require muon to stop in the TPC
///
/// \author  dcaratelli@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"
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
 
  class FilterStoppingMuon : public art::EDFilter 
  {  
  // explicit EDFilter(ParameterSet const&)  
  public:

    explicit FilterStoppingMuon(fhicl::ParameterSet const &pset);
    
    bool filter(art::Event&) ;
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
  }

  //-----------------------------------------------------------------------
  bool FilterStoppingMuon::filter(art::Event& evt) 
  {
    // for c2: interactionDesired is unused
    //bool interactionDesired(false);
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
