////////////////////////////////////////////////////////////////////////////
//
// ParticleInventory.cc
// Author: JStock
// EMail:  jason.stock@mines.sdsmt.edu
// 2017-09-12
//
////////////////////////////////////////////////////////////////////////////

//STL includes
#include <map>
//ROOT includes
//Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
//LArSoft includes
#include "larsim/MCCheater/ParticleInventory.h"
#include "nutools/ParticleNavigation/EmEveIdCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/Simulation/SimListUtils.h"
#include "lardataobj/Simulation/sim.h"
#include "lardata/Utilities/AssociationUtil.h"


namespace cheat{

  ParticleInventory::ParticleInventory( )
  {
        fG4ModuleLabel = "largeant"; //This should be replaced to make it fhicl configurable
  }

  //----------------------------------------------------------------------
  ParticleInventory::~ParticleInventory()
  {
  } 

  //-----------------------------------------------------------------------
  void ParticleInventory::ClearEvent(){
    fParticleList.clear();
    fMCTruthList.clear();
    fTrackIdToMCTruthIndex.clear();
  }

  //deliverables

  //-----------------------------------------------------------------------
  //TrackIdToParticlePtr
  const simb::MCParticle* ParticleInventory::TrackIdToParticleP(int const& id) const {
    if(!this->CanRun()){throw;}
    sim::ParticleList::const_iterator part_it = fParticleList.find(id);
    if(part_it == fParticleList.end()){
      mf::LogWarning("ParticleInventory") << "Particle with TrackId: " 
        << id << " not found in inventory. "
        << "Returning null pointer.";
      return 0;
    }
    return part_it->second;
  }//End TrackIdToParticle


  //-----------------------------------------------------------------------
  const simb::MCParticle* ParticleInventory::TrackIdToMotherParticleP(int const& id) const
  {   
    if(!this->CanRun()){throw;}
    return this->TrackIdToParticleP(fParticleList.EveId(abs(id)));
  }

  //-----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& ParticleInventory::TrackIdToMCTruthP(int const& id) const
  {
    // find the entry in the MCTruth collection for this track id
    if(!this->CanRun()){throw;}
    auto mctItr = fTrackIdToMCTruthIndex.find(abs(id));
    if(mctItr!=fTrackIdToMCTruthIndex.end()){ 
      int partIndex = mctItr->second;
      return fMCTruthList.at(partIndex);
    }else{
      throw cet::exception("ParticleInventory") << "Attempt to find MCTruth for TrackId: "
        << id <<" has failed.";
    }
  }

  //-----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& ParticleInventory::ParticleToMCTruthP(const simb::MCParticle* p) const
  {
    if(!this->CanRun()){throw;}
    return this->TrackIdToMCTruthP(p->TrackId());
  }

  //-----------------------------------------------------------------------
  const std::vector< art::Ptr<simb::MCTruth> >& ParticleInventory::MCTruthVector() const {
    if(!this->CanRun()){throw;}
    return fMCTruthList;
  }

  //-----------------------------------------------------------------------
  const std::vector<const simb::MCParticle*> ParticleInventory::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
  {
    if(!this->CanRun()){throw;}
    std::vector<const simb::MCParticle*> ret;
    // sim::ParticleList::value_type is a pair (track Id, particle pointer)
    for (const sim::ParticleList::value_type& TrackIdpair: fParticleList) {
      if( this->TrackIdToMCTruthP(TrackIdpair.first) == mct )
        ret.push_back(TrackIdpair.second);
    }
    return ret;
  }

  //-----------------------------------------------------------------------
  std::set<int> ParticleInventory::GetSetOfTrackIds() const{
    if(!this->CanRun()){throw;}
    std::set<int> ret;
    for( auto partItr=fParticleList.begin(); partItr!=fParticleList.end(); ++partItr){
      ret.emplace((partItr->second)->TrackId());
    }
    return ret;
  }

  //-----------------------------------------------------------------------
  std::set<int> ParticleInventory::GetSetOfEveIds() const{
    if(!this->CanRun()){throw;}
    std::set<int> ret;
    std::set<int> tIds=this->GetSetOfTrackIds();
    for(auto tId : tIds){
      ret.emplace(fParticleList.EveId(tId));
    }
    return ret;
  }


} //namespace



