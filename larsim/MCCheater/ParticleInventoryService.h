////////////////////////////////////////////////////////////////////////
// \file ParticleInventoryService.h
// \a service for managing the ParticleInventory when run in art.
//
// \author jason.stock@mines.sdsmt.edu
// \Based on the original BackTracker by Brian Rebel (brebel@fnal.gov)
////////////////////////////////////////////////////////////////////////
#ifndef CHEAT_PARTICLEINVENTORYSERVICESERVICE_H
#define CHEAT_PARTICLEINVENTORYSERVICESERVICE_H

#include <vector>

#include "canvas/Persistency/Common/Assns.h"


#include "larsim/MCCheater/ParticleInventory.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h" //Needed for Legacy support


#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/ParticleNavigation/EveIdCalculator.h"
#include "nusimdata/SimulationBase/MCTruth.h"



namespace cheat{
  class ParticleInventoryService
  {
    public:

      struct fhiclConfig{
        fhicl::Table<ParticleInventory::fhiclConfig> ParticleInventoryTable{
          fhicl::Name("providerConfigParticleInventory"),
          fhicl::Comment("This is the fhicl configuration for the ParticleInventory Service Provider") };
      };

      //attempting to be compliant with ServiceUtil.h. Should ask LArSoft expert to review.
      using provider_type = ParticleInventory;
      provider_type const* provider() const
      { return static_cast<provider_type const*>(&fPartInv); }
//      const ParticleInventory* provider() const {return &fPartInv;}
//      const std::shared_ptr<cheat::ParticleInventory> AccessInventory(){ return fPartInv; }

//      provider_type const* provider() const {return static_cast<provider_type const*>(fPartInv);}

      ParticleInventoryService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ParticleInventoryService(const fhiclConfig& config, art::ActivityRegistry& reg);
      ~ParticleInventoryService();

      //Move this function into the ParticleInventory.cpp file, and give it an appropriate CheckReady and Prep before the return.
      const sim::ParticleList& ParticleList() const { return fPartInv.ParticleList(); } //This should be replaced with a public struct so we can get away from the nutools dependency.
      void SetEveIdCalculator(sim::EveIdCalculator *ec) { fPartInv.SetEveIdCalculator(ec); }

      //Does this make sense? A track Id to a single particle? This is not a one to one relationship.
      const simb::MCParticle* TrackIdToParticle_P(int const& id);
      simb::MCParticle        TrackIdToParticle(int const& id)
      { return *(this->TrackIdToParticle_P(id)); }//Users are encouraged to use TrackIdToParticleP

      const simb::MCParticle* TrackIdToMotherParticle_P(int const& id);
      simb::MCParticle        TrackIdToMotherParticle(int const& id)//Users are encouraged to use TrackIdToMotherParticleP
      { return *(this->TrackIdToMotherParticle_P(id)); }

      const art::Ptr<simb::MCTruth>& TrackIdToMCTruth_P(int const& id);
      simb::MCTruth                  TrackIdToMCTruth (int const& id)//Users are encouraged to use TrackIdToMCTruthP
      { return *(this->TrackIdToMCTruth_P(id)); }

      const art::Ptr<simb::MCTruth>& ParticleToMCTruth_P(const simb::MCParticle* p); //Users are encouraged to use ParticleToMCTruthP
      simb::MCTruth                  ParticleToMCTruth (const simb::MCParticle* p)
      { return *(this->ParticleToMCTruth_P(p)); } 

      const std::vector< art::Ptr<simb::MCTruth> >& MCTruthVector_Ps() ; //I don't want this to be able to return a vector of copies. Too much chance of significant memory usage.

      const std::vector<const simb::MCParticle*> MCTruthToParticles_Ps(art::Ptr<simb::MCTruth> const& mct) ; //I don't want this to be able to return a vector of copies. Too much chance of significant memory usage.

      std::set<int> GetSetOfTrackIds();
      std::set<int> GetSetOfEveIds();



    private:

      cheat::ParticleInventory fPartInv;

      const art::Event* fEvt;

      void priv_PrepEvent        ( const art::Event& evt );
      void priv_PrepParticleList            ( );
      void priv_PrepMCTruthList             ( );
      void priv_PrepTrackIdToMCTruthIndex   ( );
      bool priv_CanRun(const art::Event& evt) const;

      bool priv_ParticleListReady()     { return  fPartInv.ParticleListReady(); }
      bool priv_MCTruthListReady()      { return  fPartInv.MCTruthListReady(); }
      bool priv_TrackIdToMCTruthReady() { return  fPartInv.TrackIdToMCTruthReady();}
  };//class ParticleInventoryService

}//namespace

DECLARE_ART_SERVICE(cheat::ParticleInventoryService, LEGACY)


#endif //CHEAT_PARTICLEINVENTORYSERVICESERVICE_H

