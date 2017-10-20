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

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Provenance/ProductToken.h"

#include "fhiclcpp/types/Atom.h"
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
      ////////////////Types/////////////////
      struct fhiclConfig{
        fhicl::Atom<art::InputTag> G4ModuleLabel{fhicl::Name("G4ModuleLabel"), fhicl::Comment("The label of the LArG4 module used to produce the art file we will be backtracking in"), "largeant"};
      };
      
      ///////////Constructor///////////////
      ParticleInventory(const fhiclConfig& config );
      ParticleInventory(const fhicl::ParameterSet& pSet );
      ~ParticleInventory();

      template<typename Evt> //Template must be decalred and defined outside of the .cpp file.
        void PrepEvent        ( const Evt& evt );

      bool ParticleListReady()     const { return !( fParticleList.empty() ); }
      bool MCTruthListReady()      const { return !( (fMCTObj.fMCTruthList).empty()  ); }
      bool TrackIdToMCTruthReady() const { return !(fMCTObj.fTrackIdToMCTruthIndex.empty());}

      template<typename Evt>    
        void PrepParticleList         (const Evt& evt ) const;
      template<typename Evt>    
        void PrepTrackIdToMCTruthIndex(const Evt& evt ) const;
      template<typename Evt>    
        void PrepMCTruthList          (const Evt& evt ) const;
      template<typename Evt> 
        void PrepMCTruthListAndTrackIdToMCTruthIndex(const Evt& evt ) const ;
      template<typename Evt>
        bool CanRun(const Evt& evt) const;

      const sim::ParticleList& ParticleList() const { return fParticleList; } 
      void SetEveIdCalculator(sim::EveIdCalculator *ec) { fParticleList.AdoptEveIdCalculator(ec); }

      const std::vector< art::Ptr<simb::MCTruth> >& MCTruthList() const { return fMCTObj.fMCTruthList;}
      
        ;
      const std::map<unsigned short, unsigned short >& TrackIdToMCTruthIndex() const { return fMCTObj.fTrackIdToMCTruthIndex; }

      void ClearEvent();

      //Does this make sense? A track Id to a single particle? This is not a one to one relationship.
      const simb::MCParticle* TrackIdToParticle_P(int const& id) const;
      simb::MCParticle        TrackIdToParticle(int const& id) const
      { return *(this->TrackIdToParticle_P(id)); }//Users are encouraged to use TrackIdToParticleP

      const simb::MCParticle* TrackIdToMotherParticle_P(int const& id) const;
      simb::MCParticle        TrackIdToMotherParticle(int const& id) const//Users are encouraged to use TrackIdToMotherParticleP
      { return *(this->TrackIdToMotherParticle_P(id)); }

      const art::Ptr<simb::MCTruth>& TrackIdToMCTruth_P(int const& id) const;
      simb::MCTruth                  TrackIdToMCTruth (int const& id) const//Users are encouraged to use TrackIdToMCTruthP
      { return *(this->TrackIdToMCTruth_P(id)); }

      //New Functions go here.
      //TrackIdToEveId.
      int TrackIdToEveTrackId(const int& tid) const { return fParticleList.EveId(tid);}

      const art::Ptr<simb::MCTruth>& ParticleToMCTruth_P(const simb::MCParticle* p) const; //Users are encouraged to use ParticleToMCTruthP
      simb::MCTruth                  ParticleToMCTruth (const simb::MCParticle* p) const
      { return *(this->ParticleToMCTruth_P(p)); } 

      const std::vector< art::Ptr<simb::MCTruth> >& MCTruthVector_Ps() const; //I don't want this to be able to return a vector of copies. Too much chance of significant memory usage.

      const std::vector<const simb::MCParticle*> MCTruthToParticles_Ps(art::Ptr<simb::MCTruth> const& mct) const; //I don't want this to be able to return a vector of copies. Too much chance of significant memory usage.

      std::set<int> GetSetOfTrackIds() const;
      std::set<int> GetSetOfEveIds() const;


    private:
      mutable sim::ParticleList                       fParticleList;
      struct MCTObjects{
        std::vector< art::Ptr<simb::MCTruth> >  fMCTruthList;   //there is some optimization that can be done here.
        std::map<unsigned short, unsigned short > fTrackIdToMCTruthIndex;
      };
      mutable MCTObjects fMCTObj;  
      //For fhicl validation, makea config struct
      art::InputTag fG4ModuleLabel;





  };//class ParticleInventory

}//namespace

#include "ParticleInventory.tpp"

#endif //CHEAT_PARTICLEINVENTORY_H

