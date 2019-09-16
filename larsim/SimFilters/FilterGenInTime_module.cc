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
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/sim.h"
#include "larcore/Geometry/Geometry.h"

// C++ Includes
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

    std::vector<std::array<double, 6>> fCryostatBoundaries; //!< boundaries of each cryostat
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
      produces< std::vector<simb::MCTruth> >("outtime"); 
    }
    std::cout << "New Filter!\n";

  }

  void FilterGenInTime::beginJob(){
    auto const& geom = *art::ServiceHandle<geo::Geometry const>();
    for (auto const &cryo: geom.IterateCryostats()) {
      std::array<double, 6> this_cryo_boundaries {};
      cryo.Boundaries(&this_cryo_boundaries[0]);
      fCryostatBoundaries.push_back(this_cryo_boundaries);
    }
  }


  bool FilterGenInTime::KeepParticle(simb::MCParticle const& part) const {
    const TLorentzVector& v4 = part.Position();
    const TLorentzVector& p4 = part.Momentum();
    // origin of particle
    double x0[3] = {v4.X(),  v4.Y(),  v4.Z() };
    // normalized direction of particle
    double dx[3] = {p4.Px() / p4.Vect().Mag(), p4.Py() / p4.Vect().Mag(), p4.Pz() / p4.Vect().Mag()};

    // tolernace for treating number as "zero"
    double eps = 1e-5;

    //check to see if particle crosses boundary of any cryostat within appropriate time window
    std::vector<bool> intersects_cryo(fCryostatBoundaries.size(), false);
    std::vector<bool> inside_cryo(fCryostatBoundaries.size(), false);
    std::vector<double> distance_to_cryo(fCryostatBoundaries.size(), 0.);

    // Check to see if particle intersects any cryostat
    //
    // Algorithmically, this is looking for ray-box intersection. This is a common problem in 
    // computer graphics. The algorithm below is taken from "Graphics Gems", Academic Press, 1990
    for (size_t i_cryo = 0; i_cryo < fCryostatBoundaries.size(); i_cryo++) {
      auto const &bound = fCryostatBoundaries[i_cryo];
      std::array<int, 3> quadrant {}; // 0 == RIGHT, 1 == LEFT, 2 == MIDDLE
      std::array<double, 3> candidatePlane {};
      std::array<double, 3> coord {};

      std::array<double, 3> bound_lo = {{bound[0], bound[2], bound[4]}};
      std::array<double, 3> bound_hi = {{bound[1], bound[3], bound[5]}};

      // First check if origin is inside box
      // Also check which of the two planes in each dimmension is the 
      // "candidate" for the ray to hit
      bool inside = true;
      for (int i = 0; i < 3; i++) {
        if (x0[i] < bound_lo[i]) {
          quadrant[i] = 1; // LEFT
          candidatePlane[i] = bound_lo[i];
          inside = false;
        }
        else if (x0[i] > bound_hi[i]) {
          quadrant[i] = 0; // RIGHT
          candidatePlane[i] = bound_hi[i];
          inside = false;
        }
        else {
          quadrant[i] = 2; // MIDDLE
        }
      }

      if (inside) {
        inside_cryo[i_cryo] = true;
        // if we're inside the cryostat, then we do intersect it
        intersects_cryo[i_cryo] = true;
        continue;
      }

      // ray origin is outside the box -- calculate the distance to the cryostat and see if it intersects
      
      // calculate distances to candidate planes
      std::array<double, 3> maxT {};
      for (int i = 0; i < 3; i++) {
        if (quadrant[i] != 2 /* MIDDLE */ && abs(dx[i]) > eps) {
          maxT[i] = (candidatePlane[i] - x0[i]) / dx[i];
        }
        // if a ray origin is between two the two planes in a dimmension, it would never hit that plane first
        else {
          maxT[i] = -1;
        }
      }

      // The plane on the box that the ray hits is the one with the largest distance
      int whichPlane = 0;
      for (int i = 1; i < 3; i++) {
        if (maxT[whichPlane] < maxT[i]) whichPlane = i;
      }

      // check if the candidate intersection point is inside the box

      // no intersection
      if (maxT[whichPlane] < 0.) {
        intersects_cryo[i_cryo] = false;
        continue;
      }

      for (int i = 0; i < 3; i++) {
        if (whichPlane != i) {
          coord[i] = x0[i] + maxT[whichPlane] * dx[i];
        }
        else {
          coord[i] = candidatePlane[i];
        }
      }


      // check if intersection is in box
      intersects_cryo[i_cryo] = true;
      for (int i = 0; i < 3; i++) {
        if (coord[i] < bound_lo[i] || coord[i] > bound_hi[i]) {
          intersects_cryo[i_cryo] = false;
        }
      }

      if (intersects_cryo[i_cryo]) {
        distance_to_cryo[i_cryo] = maxT[whichPlane];
      }
    }

    // check if any cryostats are intersected in-time
    for (size_t i_cryo = 0; i_cryo < fCryostatBoundaries.size(); i_cryo++) {

      // If the particle originates inside the cryostat, then
      // we can't really say when it will leave. Thus, accept
      // the particle
      if (inside_cryo[i_cryo]) {
        return true;
      }
      // otherwise check arrival time at boundary of cryostat
      if (intersects_cryo[i_cryo]){
        double ptime = (distance_to_cryo[i_cryo] * 1e-2 /* cm -> m */) / (TMath::C()*sqrt(1-pow(part.Mass()/part.E(),2))) /* velocity */;
        double totT=part.T()+ptime*1e9 /* s -> ns */;
        if(totT>fMinT && totT<fMaxT){
          return true;
        }
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
    art::ServiceHandle<geo::Geometry const> geom;

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
