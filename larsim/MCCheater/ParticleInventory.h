////////////////////////////////////////////////////////////////////////
// ParticleInventory.h
// Provide a single interface for building and accessing truth information from events for backtracking services.
//
// author jason.stock@mines.sdsmt.edu
// Based on the original BackTracker by Brian Rebel (brebel@fnal.gov)
//
// This module may look strange at first glance beacause of the mutable 
// object in const functions whose sole purpose is to change the mutable 
// objects.
// This is done because the returns from the ParticleInventory really 
// should be const, and anything called from them must then also be const,
// and finally we get to the mutables. These are cached objects to prevent 
// repeated and costly access to objects in the event.
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//DOXYGEN DOCUMENTATION
////////////////////////////////////////////////////////////////////////
/**
 * \file ParticleInventory.h
 * \author jason.stock@mines.sdsmt.edu  
 * \brief Header for the ParticleInvenotry Service Provider.
 *
 * The ParticleInventory is an art independent service provider 
 * for retreiving truth information about tracks and particles.
 * The ParticleInventory, BackTracker, and PhotonBackTracker make 
 * a complete toolset for retreiving truth information from an event.
 */

/** \struct cheat::ParticleInventory::ParticleInventoryConfig
 *  \brief FHICL Validation Object
 *  This struct is used for loading the fhicl configuration.
 */ 
/** \var cheat::ParticleInventory::ParticleInventoryConfig::G4ModuleLabel
 *  \brief An atom.
 *  FHICL Atom for retreiving the module label to be used in retreiving
 *  information from the art event.
 *  */
/** \struct cheat::ParticleInventory::MCTObjects
 *  \brief  A simple struct to contain the MC Truth information.
 *  \var cheat::ParticleInventory::MCTObjects::fMCTruthList;
 *  \brief A vector containing the MCTruth objects
 *  \var cheat::ParticleInventory::MCTobjects::fTrackIdToMCTruthIndex
 *  \brief a map linking trackIds to the location of MCTruth for fast lookup.
 *  */

/** \fn void cheat::ParticleInventory::PrepEvent        ( const Evt& evt )
 *  \brief Function to set up the ParticleInventory state for an event.
 *  This is a function to tell the ParticleInventory to prepare itself to work with a particular event.
 *  @param evt
 *  \brief The event the ParticleInventory should work with. *Note. This use breaks the multithreading model because the service has a "state".
 *  */
 
/** \fn bool cheat::ParticleInventory::ParticleListReady()     const
 *  \brief A simple check to determine if the ParticleList has already been prepared for this event or not.
 *  */
/** \fn bool cheat::ParticleInventory::MCTruthListReady()      const 
 *  \brief A simple check to determine if the MCTruthList has already been prepared and cached or not.
 *  */
/** \fn bool cheat::ParticleInventory::TrackIdToMCTruthReady() const
 *  \brief A simple check to determine if the TrackIdToMCTruth map has been prepared or not.
 *  */
/** \fn void cheat::ParticleInventory::PrepParticleList         (const Evt& evt ) const
 *  \brief A function to load the ParticleList and cache it
 *  This function will find the particle list and load it for later use. 
 *  Ideally this would would be used for a "lazy" loading of the backtracker, 
 *  but this does not work in the current setup of art.
 *  */
/** \fn void cheat::ParticleInventory::PrepTrackIdToMCTruthIndex(const Evt& evt ) const
 *  \breif A function to prepare and cache a map of TrackIds and MCTruth object indicies from fMCTruthList.
 *  */
/** \fn void cheat::ParticleInventory::PrepMCTruthList          (const Evt& evt ) const
 *  \brief A function to load and cache the MCTruthList of the event.
 *  */
/** \fn void cheat::ParticleInventory::PrepMCTruthListAndTrackIdToMCTruthIndex(const Evt& evt ) const
 *  \brief A function to make both PrepTrackIdToMCTruthIndex and PrepMCTruthList run when both are needed.
 *  */
/** \fn bool cheat::ParticleInventory::CanRun(const Evt& evt) const
 *  \brief A short function to check if use of the backtracker is appropriate or not based on the type of input file.
 *  This function simply checks to see if the file loaded is real data, or MC Simulation, as backtracking on real data makes no sense.
 *  If one does try to backtrack real data, this will throw and exception.
 *  */
/** \fn const sim::ParticleList& ParticleList() const
 *  \breif Get the ParticleList from an event.
 *  */
/** \fn  void SetEveIdCalculator(sim::EveIdCalculator *ec)
 *  \brief Set the EveIdCalculator to use for this ParticleList. 
 *  Set the EveIdCalculator to use for this ParticleList. If you are 
 *  going to over-ride the default EveIdCalculator, you must do it for 
 *  EVERY event (as a new particle list is adopted for each event.
 */
#ifndef CHEAT_PARTICLEINVENTORY_H
#define CHEAT_PARTICLEINVENTORY_H

#include <vector>

#include "art/Framework/Principal/Handle.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindOneP.h"
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
      struct ParticleInventoryConfig{
        fhicl::Atom<art::InputTag> G4ModuleLabel{
          fhicl::Name("G4ModuleLabel"), 
          fhicl::Comment("The label of the LArG4 module used to produce the art file we will be examining"), 
          "largeant"};
      };

      //using provider_type = ParticleInventory;
      //cheat::ParticleInventory const* provider() const
      //{ return static_cast<cheat::ParticleInventory const*>(this); }
      
      ///////////Constructor///////////////
      ParticleInventory(const ParticleInventoryConfig& config );
      ParticleInventory(const fhicl::ParameterSet& pSet );
      ParticleInventory(ParticleInventory const&) = delete;

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
      const std::map< int,  int >& TrackIdToMCTruthIndex() const { return fMCTObj.fTrackIdToMCTruthIndex; }

      void ClearEvent();

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
        std::map< int,  int > fTrackIdToMCTruthIndex;
      };
      mutable MCTObjects fMCTObj;  
      //For fhicl validation, makea config struct
      art::InputTag fG4ModuleLabel;





  };//class ParticleInventory

}//namespace

#include "ParticleInventory.tcc"

#endif //CHEAT_PARTICLEINVENTORY_H

