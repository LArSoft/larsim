////////////////////////////////////////////////////////////////////////
/// \file  FilterGenInTime_module.cc
/// \brief EDFilter to require projected generator trajectories in volumes within a particular time window.
///
/// \author  Matthew.Bass@physics.ox.ac.uk
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

// Root includes
#include <TMath.h>

namespace simfilter {  
 
  class FilterGenInTime : public art::EDFilter 
  {  
  public:

    explicit FilterGenInTime(fhicl::ParameterSet const &pset);
    
    bool filter(art::Event&) ;
      
    virtual void beginJob();
    
    
  private:
    bool KeepParticle(simb::MCParticle const& part) const;
    
    double fbounds[6] = {0.};
    double fMinKE; //<only keep based on particles with greater than this energy   
    bool   fKeepOnlyMuons; //keep based only on muons if enabled
    double fMinT,fMaxT; //<time range in which to keep particles
    bool   fSortParticles; //create new MCTruth collections with particles sorted by their timing
    bool   fAlwaysPass; // flag to have filter always pass (to be used with sorting...)
  };

} // namespace simfilter

namespace simfilter {

  FilterGenInTime::FilterGenInTime(fhicl::ParameterSet const& pset) :
    EDFilter{pset},
    fMinKE    (pset.get< double > ("MinEnergy"  , 0.0)       )
    , fKeepOnlyMuons    (pset.get< bool > ("KeepOnlyMuons", false)      )
    , fMinT    (pset.get< double > ("MinT",0.0)      )
    , fMaxT    (pset.get< double > ("MaxT")      )
    , fSortParticles ( pset.get< bool > ("SortParticles",false) )
    , fAlwaysPass (pset.get<bool>("AlwaysPass",false))
  {
    if(fSortParticles) {
      produces< std::vector<simb::MCTruth> >("intime");
      produces< std::vector<simb::MCTruth> >("outtime"); }
  }

  void FilterGenInTime::beginJob(){
    auto const& geom = *art::ServiceHandle<geo::Geometry>();
	  
	  geom.CryostatBoundaries(fbounds, 0);
  }
  
    
  bool FilterGenInTime::KeepParticle(simb::MCParticle const& part) const {
    const TLorentzVector& v4 = part.Position();
    const TLorentzVector& p4 = part.Momentum();
    double x0[3] = {v4.X(),  v4.Y(),  v4.Z() };
    double dx[3] = {p4.Px(), p4.Py(), p4.Pz()};
    
    //check to see if particle crosses boundary of cryostat within appropriate time window
    bool intersects_cryo = false;
    for (int bnd=0; bnd!=6; ++bnd) {
      if (bnd<2) {
	double p2[3] = {fbounds[bnd],  x0[1] + (dx[1]/dx[0])*(fbounds[bnd] - x0[0]), x0[2] + (dx[2]/dx[0])*(fbounds[bnd] - x0[0])};
	if ( p2[1] >= fbounds[2] && p2[1] <= fbounds[3] && 
	     p2[2] >= fbounds[4] && p2[2] <= fbounds[5] ) {
	  intersects_cryo = true;
	  break;
	}
      }
      else if (bnd>=2 && bnd<4) {
	double p2[3] = {x0[0] + (dx[0]/dx[1])*(fbounds[bnd] - x0[1]), fbounds[bnd], x0[2] + (dx[2]/dx[1])*(fbounds[bnd] - x0[1])};
	if ( p2[0] >= fbounds[0] && p2[0] <= fbounds[1] && 
	     p2[2] >= fbounds[4] && p2[2] <= fbounds[5] ) {
	  intersects_cryo = true;
	  break;
	}
      }
      else if (bnd>=4) {
	double p2[3] = {x0[0] + (dx[0]/dx[2])*(fbounds[bnd] - x0[2]), x0[1] + (dx[1]/dx[2])*(fbounds[bnd] - x0[2]), fbounds[bnd]};
	if ( p2[0] >= fbounds[0] && p2[0] <= fbounds[1] && 
	     p2[1] >= fbounds[2] && p2[1] <= fbounds[3] ) {
	  intersects_cryo = true;
	  break;
	}
      }
    }
    
    if (intersects_cryo){
      //check its arrival time at the origin(base time + propagation time)
      double d=sqrt((x0[0]*x0[0] + x0[1]*x0[1] + x0[2]*x0[2]));
      double ptime=d/(100.*TMath::C()*sqrt(1-pow(part.Mass()/part.E(),2)));
      double totT=part.T()+ptime*1e9;
      if(totT>fMinT && totT<fMaxT){
        return true;
      }
    }
    return false;
  }
  
  bool FilterGenInTime::filter(art::Event& evt){

    std::unique_ptr< std::vector<simb::MCTruth> > truthInTimePtr(new std::vector<simb::MCTruth>(1));
    std::unique_ptr< std::vector<simb::MCTruth> > truthOutOfTimePtr(new std::vector<simb::MCTruth>(1));

    simb::MCTruth & truthInTime = truthInTimePtr->at(0);
    simb::MCTruth & truthOutOfTime = truthOutOfTimePtr->at(0);

    
    //get the list of particles from this event
    art::ServiceHandle<geo::Geometry> geom;
    
    std::vector< art::Handle< std::vector<simb::MCTruth> > > allmclists;
    evt.getManyByType(allmclists);

    bool keepEvent=false;
    for(size_t mcl = 0; mcl < allmclists.size(); ++mcl){
      art::Handle< std::vector<simb::MCTruth> > mclistHandle = allmclists[mcl];
      for(size_t m = 0; m < mclistHandle->size(); ++m){
	art::Ptr<simb::MCTruth> mct(mclistHandle, m);
	
	for(int ipart=0;ipart<mct->NParticles();ipart++){
	  bool kp=KeepParticle(mct->GetParticle(ipart));

	  if(kp
	     && (!fKeepOnlyMuons || abs(mct->GetParticle(ipart).PdgCode())==13 )
	     && mct->GetParticle(ipart).E()-mct->GetParticle(ipart).Mass()>fMinKE){
	    keepEvent = true;
	    if(!fSortParticles) break;
	  }

	  if(fSortParticles){
	    simb::MCParticle particle = mct->GetParticle(ipart);
	    if(kp) truthInTime.Add(particle);
	    if(!kp) truthOutOfTime.Add(particle);
	  }
	  
	}//end loop over particles

	if(!fSortParticles && keepEvent) break;

      }//end loop over mctruth col

      if(!fSortParticles && keepEvent) break;
      
    }//end loop over all mctruth lists

    if(fSortParticles){
      evt.put(std::move(truthInTimePtr),"intime");
      evt.put(std::move(truthOutOfTimePtr),"outtime");
    }
    
    return (keepEvent || fAlwaysPass);
  }
  
} // namespace simfilter

namespace simfilter {
  
  DEFINE_ART_MODULE(FilterGenInTime)
  
} // namespace simfilter
