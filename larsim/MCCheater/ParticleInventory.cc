////////////////////////////////////////////////////////////////////////
//
// ParticleInventory.cc
// Author: JStock
// EMail:  jason.stock@mines.sdsmt.edu
// 2017-09-12
//
////////////////////////////////////////////////////////////////////////

//STL includes
//ROOT includes
//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//LArSoft includes
#include "larsim/MCCheater/ParticleInventory.h"
#include "nusimdata/SimulationBase/MCParticle.h"


namespace cheat{

  ParticleInventory::ParticleInventory(const ParticleInventoryConfig& config )
    :fG4ModuleLabel(config.G4ModuleLabel()),
    fEveIdCalculator(config.EveIdCalculator()),
    fOverrideRealData(config.OverrideRealData())
  {
  }

  //----------------------------------------------------------------------
  ParticleInventory::ParticleInventory(const fhicl::ParameterSet& pSet )
    :fG4ModuleLabel(pSet.get<art::InputTag>("G4ModuleLabel", "largeant")),
    fEveIdCalculator(pSet.get<std::string>("EveIdCalculator", "EmEveIdCalculator")),
    fOverrideRealData(pSet.get<bool>("OverrideRealData", false))
  {
  }

  //-----------------------------------------------------------------------
  void ParticleInventory::ClearEvent(){
    fParticleList.clear();
    fMCTObj.fMCTruthList.clear();
    fMCTObj.fTrackIdToMCTruthIndex.clear();
  }

  //deliverables

  //-----------------------------------------------------------------------
  //TrackIdToParticlePtr
  const simb::MCParticle* ParticleInventory::TrackIdToParticle_P(int const& id) const {
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
  const simb::MCParticle* ParticleInventory::TrackIdToMotherParticle_P(int const& id) const
  {
    return this->TrackIdToParticle_P(fParticleList.EveId(abs(id)));
  }

  //-----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& ParticleInventory::TrackIdToMCTruth_P(int const& id) const
  {
    // find the entry in the MCTruth collection for this track id
    auto mctItr = fMCTObj.fTrackIdToMCTruthIndex.find(abs(id));
    if(mctItr!=fMCTObj.fTrackIdToMCTruthIndex.end()){
      int partIndex = mctItr->second;
      return fMCTObj.fMCTruthList.at(partIndex);
    }else{
      throw cet::exception("ParticleInventory") << "Attempt to find MCTruth for TrackId: "
        << id <<" has failed.";
    }
  }

  //-----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& ParticleInventory::ParticleToMCTruth_P(const simb::MCParticle* p) const
  {
    return this->TrackIdToMCTruth_P(p->TrackId());
  }

  //-----------------------------------------------------------------------
  const std::vector< art::Ptr<simb::MCTruth> >& ParticleInventory::MCTruthVector_Ps() const {
    return fMCTObj.fMCTruthList;
  }

  //-----------------------------------------------------------------------
  std::vector<const simb::MCParticle*> ParticleInventory::MCTruthToParticles_Ps(art::Ptr<simb::MCTruth> const& mct) const
  {
    std::vector<const simb::MCParticle*> ret;
    // sim::ParticleList::value_type is a pair (track Id, particle pointer)
    for (const sim::ParticleList::value_type& TrackIdpair: fParticleList) {
      if( this->TrackIdToMCTruth_P(TrackIdpair.first) == mct )
        ret.push_back(TrackIdpair.second);
    }
    return ret;
  }

  //-----------------------------------------------------------------------
  std::set<int> ParticleInventory::GetSetOfTrackIds() const{
    std::set<int> ret;
    for( auto partItr=fParticleList.begin(); partItr!=fParticleList.end(); ++partItr){
      ret.emplace((partItr->second)->TrackId());
    }
    return ret;
  }

  //-----------------------------------------------------------------------
  std::set<int> ParticleInventory::GetSetOfEveIds() const{
    std::set<int> ret;
    std::set<int> tIds=this->GetSetOfTrackIds();
    for(auto tId : tIds){
      ret.emplace(fParticleList.EveId(tId));
    }
    return ret;
  }


} //namespace
