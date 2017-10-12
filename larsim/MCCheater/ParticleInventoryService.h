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
      ParticleInventoryService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~ParticleInventoryService();

      const sim::ParticleList& ParticleList() const { return fPartInv.ParticleList(); } //This should be replaced with a public struct so we can get away from the nutools dependency.
      void SetEveIdCalculator(sim::EveIdCalculator *ec) { fPartInv.SetEveIdCalculator(ec); }

      //Does this make sense? A track Id to a single particle? This is not a one to one relationship.
      const simb::MCParticle* TrackIdToParticleP(int const& id);
      simb::MCParticle        TrackIdToParticle(int const& id)
      { return *(this->TrackIdToParticleP(id)); }//Users are encouraged to use TrackIdToParticleP

      const simb::MCParticle* TrackIdToMotherParticleP(int const& id);
      simb::MCParticle        TrackIdToMotherParticle(int const& id)//Users are encouraged to use TrackIdToMotherParticleP
      { return *(this->TrackIdToMotherParticleP(id)); }

      const art::Ptr<simb::MCTruth>& TrackIdToMCTruthP(int const& id);
      simb::MCTruth                  TrackIdToMCTruth (int const& id)//Users are encouraged to use TrackIdToMCTruthP
      { return *(this->TrackIdToMCTruthP(id)); }

      const art::Ptr<simb::MCTruth>& ParticleToMCTruthP(const simb::MCParticle* p); //Users are encouraged to use ParticleToMCTruthP
      simb::MCTruth                  ParticleToMCTruth (const simb::MCParticle* p)
      { return *(this->ParticleToMCTruthP(p)); } 

      const std::vector< art::Ptr<simb::MCTruth> >& MCTruthVector() ; //I don't want this to be able to return a vector of copies. Too much chance of significant memory usage.

      const std::vector<const simb::MCParticle*> MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) ; //I don't want this to be able to return a vector of copies. Too much chance of significant memory usage.

      std::set<int> GetSetOfTrackIds();
      std::set<int> GetSetOfEveIds();


    private:
      std::string fG4ModuleLabel;
      cheat::ParticleInventory  fPartInv;

      const art::Event* fEvt;

      void priv_PrepEvent        ( const art::Event& evt );
      void priv_PrepParticleList            ( );
      void priv_PrepMCTruthList             ( );
      void priv_PrepTrackIdToMCTruthIndex   ( );
      bool priv_CanRun();

      bool priv_ParticleListReady()     { return  fPartInv.ParticleListReady(); }
      bool priv_MCTruthListReady()      { return  fPartInv.MCTruthListReady(); }
      bool priv_TrackIdToMCTruthReady() { return  fPartInv.TrackIdToMCTruthReady();}
  };//class ParticleInventoryService

}//namespace


#endif //CHEAT_PARTICLEINVENTORYSERVICESERVICE_H

