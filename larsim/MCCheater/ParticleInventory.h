////////////////////////////////////////////////////////////////////////
// \file ParticleInventory.h
// \brief Provide a single interface for building and accessing truth information from events for backtracking services.
//
// \author jason.stock@mines.sdsmt.edu
// \Based on the original BackTracker by Brian Rebel (brebel@fnal.gov)
//
// This module may look strange at first glance beacause of the mutable object in const functions whose sole purpose is to change the mutable objects.
// This is done because the returns from the ParticleInventory really should be const, and anything called from them must then also be const,
// and finally we get to the mutables. These are cached objects to prevent repeated and costly access to objects in the event.
////////////////////////////////////////////////////////////////////////
#ifndef CHEAT_PARTICLEINVENTORY_H
#define CHEAT_PARTICLEINVENTORY_H

#include <vector>

#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/Assns.h"

#include "fhiclcpp/ParameterSet.h"
#include "nutools/ParticleNavigation/EmEveIdCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/ParticleNavigation/EveIdCalculator.h"
#include "nusimdata/SimulationBase/MCTruth.h"



namespace cheat{
  class ParticleInventory
  {
    public:
      ParticleInventory();
      ~ParticleInventory();

      template<typename Evt> //Template must be decalred and defined outside of the .cpp file.
        void PrepEvent        ( const Evt& evt );

      bool CanRun()                const { return fCanRun;}
      bool ParticleListReady()     const { return !( fParticleList.empty() ); }
      bool MCTruthListReady()      const { return !( fMCTruthList.empty()  ); }
      bool TrackIdToMCTruthReady() const { return !(fTrackIdToMCTruthIndex.empty());}

      template<typename Evt>    
        void PrepParticleList(const Evt& evt ) const;
      template<typename Evt>    
        void PrepTrackIdToMCTruthIndex(Evt evt ) const;
      template<typename Evt>    
        void PrepMCTruthList             (Evt& evt ) const;
      template<typename Evt>    
        void CheckCanRun( Evt& evt);

      const sim::ParticleList& ParticleList() const { return fParticleList; } //This should be replaced with a public struct so we can get away from the nutools dependency.
      void SetEveIdCalculator(sim::EveIdCalculator *ec) { fParticleList.AdoptEveIdCalculator(ec); }

      void ClearEvent();

      //Does this make sense? A track Id to a single particle? This is not a one to one relationship.
      const simb::MCParticle* TrackIdToParticleP(int const& id) const;
      simb::MCParticle        TrackIdToParticle(int const& id) const
      { return *(this->TrackIdToParticleP(id)); }//Users are encouraged to use TrackIdToParticleP

      const simb::MCParticle* TrackIdToMotherParticleP(int const& id) const;
      simb::MCParticle        TrackIdToMotherParticle(int const& id) const//Users are encouraged to use TrackIdToMotherParticleP
      { return *(this->TrackIdToMotherParticleP(id)); }

      const art::Ptr<simb::MCTruth>& TrackIdToMCTruthP(int const& id) const;
      simb::MCTruth                  TrackIdToMCTruth (int const& id) const//Users are encouraged to use TrackIdToMCTruthP
      { return *(this->TrackIdToMCTruthP(id)); }

      const art::Ptr<simb::MCTruth>& ParticleToMCTruthP(const simb::MCParticle* p) const; //Users are encouraged to use ParticleToMCTruthP
      simb::MCTruth                  ParticleToMCTruth (const simb::MCParticle* p) const
      { return *(this->ParticleToMCTruthP(p)); } 

      const std::vector< art::Ptr<simb::MCTruth> >& MCTruthVector() const; //I don't want this to be able to return a vector of copies. Too much chance of significant memory usage.

      const std::vector<const simb::MCParticle*> MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const; //I don't want this to be able to return a vector of copies. Too much chance of significant memory usage.

      std::set<int> GetSetOfTrackIds() const;
      std::set<int> GetSetOfEveIds() const;


    private:
      bool                                    fCanRun=0;
      mutable sim::ParticleList                       fParticleList;
      mutable std::vector< art::Ptr<simb::MCTruth> >  fMCTruthList;   //there is some optimization that can be done here.
      mutable std::map<unsigned short, unsigned short > fTrackIdToMCTruthIndex;
      std::string fG4ModuleLabel;




  };//class ParticleInventory

}//namespace

#include "ParticleInventory.tpp"

#endif //CHEAT_PARTICLEINVENTORY_H

