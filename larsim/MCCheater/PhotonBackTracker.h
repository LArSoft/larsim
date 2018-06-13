////////////////////////////////////////////////////////////////////
//
// \file PhotonBackTracker.h
// \brief The functions needed for the PhotonBackTracker class needed by the PhotonBackTrackerService in order to connect truth information with reconstruction.
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
///////////////////////////////////////////////////////////////////
//
//TODO: Impliment alternate backtracking scheme developed by T. Usher
//
///////////////////////////////////////////////////////////////////
#ifndef CHEAT_PHOTONBACKTRACKER_H
#define CHEAT_PHOTONBACKTRACKER_H

//CPP
#include <map>
#include <vector>

//Framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

//LArSoft
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/ProviderPack.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/MCCheater/ParticleInventory.h"



namespace cheat{
  class PhotonBackTracker
  {
    public:

      struct fhiclConfig{
        fhicl::Atom<double> Delay{fhicl::Name("Delay"), fhicl::Comment("The delay time needed to correctly account for the optical simulation and optical systems simulation. (The time between when a g4partcile was made, and when the simulation write out says a signal was recorded)."), 0};
        fhicl::Atom<art::InputTag> G4ModuleLabel{fhicl::Name("G4ModuleLabel"), fhicl::Comment("The label of the LArG4 module used to produce the art file we will be using."), "largeant"};
        fhicl::Atom<art::InputTag> OpHitLabel{fhicl::Name("OpHitLabel"), fhicl::Comment("The default label for the module to use when grabbing OpHits"), "ophit"}; //This should be removed and replaced with some way to access the OpHitLabel given by the user in their own analysis module to avoid differing definitions.
        fhicl::Atom<art::InputTag> OpFlashLabel{fhicl::Name("OpFlashLabel"), fhicl::Comment("The default label for the module to use when grabbing OpFlash"), "opflash"}; //This should be removed and replaced with some way to access the OpFlashLabel given by the user in their own analysis module to avoid differing definitions.
        fhicl::Atom<double> MinOpHitEnergyFraction{fhicl::Name("MinOpHitEnergyFraction"), fhicl::Comment("The minimum contribution an energy deposit must make to a Hit to be considered part of that hit."),0.010};
      };

      //----------------------------------------------------------------
      PhotonBackTracker(fhiclConfig const& config,
          const cheat::ParticleInventory* partInv,
          const geo::GeometryCore* geom);//,
//          const detinfo::DetectorClocks* detClock);

      //----------------------------------------------------------------
      PhotonBackTracker( fhicl::ParameterSet const& pSet,
          const cheat::ParticleInventory* partInv,
          const geo::GeometryCore* geom);//,
//          const detinfo::DetectorClocks* detClock);

      //----------------------------------------------------------------
      PhotonBackTracker(PhotonBackTracker const&) = delete;

      //-----------------------------------------------------
      template<typename Evt>
        void PrepEvent( Evt const& evt);

      //-----------------------------------------------------
      template<typename Evt>
        const bool CanRun( Evt const& evt) ;

      //-----------------------------------------------------
      template<typename Evt>
        void PrepOpDetBTRs(Evt const& evt);

      //-----------------------------------------------------
      template<typename Evt>
        void PrepOpFlashToOpHits(Evt const& evt);

      //----------------------------------------------------- /*NEW*/
      const std::vector<art::Ptr<recob::OpHit>> OpFlashToOpHits_Ps(art::Ptr<recob::OpFlash>& flash_P) const;

      //----------------------------------------------------- /*NEW*/
      const std::vector<double> OpFlashToXYZ(art::Ptr<recob::OpFlash>& flash_P) const ;

      //----------------------------------------------------- /*NEW*/
      const std::set<int> OpFlashToTrackIds(art::Ptr<recob::OpFlash>& flash_P) const;

      //-----------------------------------------------------
      void ClearEvent();

      //-----------------------------------------------------
      const bool BTRsReady() ;

      //-----------------------------------------------------
      const bool OpFlashToOpHitsReady();

      //-----------------------------------------------------
      std::vector<art::Ptr<sim::OpDetBacktrackerRecord>> const& OpDetBTRs() ;

      //-----------------------------------------------------
      const std::vector< const sim::SDP* > TrackIdToSimSDPs_Ps(int const& id) ;

      //-----------------------------------------------------
      const std::vector< const sim::SDP* > TrackIdToSimSDPs_Ps(int const& id,  geo::View_t const& view) ;

      //-----------------------------------------------------
      const art::Ptr< sim::OpDetBacktrackerRecord > FindOpDetBTR(int const& opDetNum) const;

      //-----------------------------------------------------
      const std::vector < sim::TrackSDP > OpDetToTrackSDPs(int const& OpDetNum, double const& opHit_start_time, double const& opHit_end_time) const;

      //-----------------------------------------------------
      const std::vector< sim::TrackSDP> OpHitToTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P) const ;

      //-----------------------------------------------------
      const std::vector< sim::TrackSDP> OpHitToTrackSDPs(recob::OpHit const& opHit) const ;

      //-----------------------------------------------------
      const std::vector < int > OpHitToTrackIds(recob::OpHit const& opHit) const;

      //-----------------------------------------------------
      const std::vector < int > OpHitToTrackIds(art::Ptr<recob::OpHit> const& opHit) const;

      //-----------------------------------------------------
      const std::vector < int > OpHitToEveTrackIds(recob::OpHit const& opHit) ;

      //-----------------------------------------------------
      const std::vector < int > OpHitToEveTrackIds(art::Ptr<recob::OpHit> const& opHit_P);

      //-----------------------------------------------------
      const std::vector< sim::TrackSDP> OpHitToEveTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P) const;

      //-----------------------------------------------------
      const std::vector< sim::TrackSDP> OpHitToEveTrackSDPs( recob::OpHit const& opHit) const ;

      //-----------------------------------------------------
      const std::vector<art::Ptr<recob::OpHit>> TrackIdToOpHits_Ps( int const& tkId, std::vector<art::Ptr<recob::OpHit>> const& hitsIn) ;

      //-----------------------------------------------------
      const std::vector<std::vector<art::Ptr<recob::OpHit>>> TrackIdsToOpHits_Ps( std::vector< int> const& tkIds, std::vector<art::Ptr<recob::OpHit>> const& hitsIn) ;

      //-----------------------------------------------------
      const std::vector< const sim::SDP* > OpHitToSimSDPs_Ps( recob::OpHit const& opHit) const ;

      //-----------------------------------------------------
      const std::vector< const sim::SDP* > OpHitToSimSDPs_Ps( art::Ptr<recob::OpHit> const& opHit_P) const;
      //
      //-----------------------------------------------------
      //      const std::vector< const sim::SDP* > OpHitsToSimSDPs_Ps( const std::vector< art::Ptr < recob::OpHit > >& opHits_Ps) ;
      const std::vector< const sim::SDP* > OpHitsToSimSDPs_Ps( std::vector< art::Ptr < recob::OpHit  > > const& opHits_Ps) const;

      //-----------------------------------------------------
      const std::vector< double > SimSDPsToXYZ(std::vector<sim::SDP> const& sdps) const&;

      //-----------------------------------------------------
      const std::vector< double > SimSDPsToXYZ(std::vector<const sim::SDP*> const& sdps_Ps ) const&;

      //-----------------------------------------------------
      const std::vector< double > OpHitToXYZ(art::Ptr<recob::OpHit> const& opHit) ;

      //-----------------------------------------------------
      const std::vector< double > OpHitToXYZ(recob::OpHit const& opHit) ;

      //-----------------------------------------------------
      const std::vector< double > OpHitsToXYZ( std::vector < art::Ptr < recob::OpHit > > const& opHits_Ps) const; /*NEW*/

      //----------------------------------------------------- /*NEW*/
      const std::unordered_set<const sim::SDP*> OpHitToEveSimSDPs_Ps(recob::OpHit const& opHit);

      //----------------------------------------------------- /*NEW*/
      const std::unordered_set<const sim::SDP*> OpHitToEveSimSDPs_Ps(art::Ptr<recob::OpHit>& opHit_P);

      //-----------------------------------------------------
      const std::set< int> GetSetOfEveIds() const ;

      //-----------------------------------------------------
      const std::set< int> GetSetOfTrackIds() const ;

      //-----------------------------------------------------
      const std::set< int> GetSetOfEveIds(std::vector< art::Ptr<recob::OpHit> > const& opHits) const ;

      //----------------------------------------------------- /*NEW*/
      const std::set< int> GetSetOfEveIds(std::vector< recob::OpHit > const& opHits) const ;

      //-----------------------------------------------------
      const std::set< int> GetSetOfTrackIds(std::vector< art::Ptr<recob::OpHit> > const& opHits) const;

      //----------------------------------------------------- /*NEW*/
      const std::set< int> GetSetOfTrackIds(std::vector< recob::OpHit > const& opHits) const;

      //-----------------------------------------------------
      const double OpHitCollectionPurity(std::set<int> const& tkIds, std::vector< art::Ptr<recob::OpHit> > const& opHits) ;

      //-----------------------------------------------------
      const double OpHitLightCollectionPurity(std::set<int> const& tkIds, std::vector< art::Ptr<recob::OpHit> > const& opHits) ;

      //-----------------------------------------------------
      const double OpHitCollectionEfficiency(std::set<int> const& tkIds,
          std::vector< art::Ptr< recob::OpHit > > const& opHits,
          std::vector< art::Ptr< recob::OpHit > > const& opHitsIn) ;

      //-----------------------------------------------------
      const double OpHitCollectionEfficiency(std::set<int> const& tkIds,
          std::vector< art::Ptr< recob::OpHit > > const& opHits,
          std::vector< art::Ptr< recob::OpHit > > const& opHitsIn,
          geo::View_t const& view)  ;

      //-----------------------------------------------------
      const double OpHitLightCollectionEfficiency(std::set<int> const& tkIds,
          std::vector< art::Ptr< recob::OpHit > > const& opHits,
          std::vector< art::Ptr< recob::OpHit > > const& opHitsIn) ;

      //-----------------------------------------------------
      const double OpHitLightCollectionEfficiency(std::set<int> const& tkIds,
          std::vector< art::Ptr< recob::OpHit > > const& opHits,
          std::vector< art::Ptr< recob::OpHit > > const& opHitsIn,
          geo::View_t const& view);

      //----------------------------------------------------- /*NEW*/
      //std::vector<sim::TrackSDP> OpFlashToTrackSDPs(art::Ptr<recob::OpFlash> flash_P);
      //----------------------------------------------------- /*NEW*/
      //std::vector<sim::TrackSDP> OpFlashToEveTrackSDPs(recob::OpFlash flash);
      //----------------------------------------------------- /*NEW*/
      //std::vector<sim::TrackSDP> OpFlashToEveTrackSDPs(art::Ptr<recob::OpFlash> flash_P);
      //----------------------------------------------------- /*NEW*/
      //std::vector<sim::SDP*> OpFlashToSimSDPs_Ps(art::Ptr<recob::OpFlash> flash_P);


    private:
      const cheat::ParticleInventory* fPartInv; //The constructor needs to put something in here
      const geo::GeometryCore* fGeom;
//      const detinfo::DetectorClocks* fDetClocks;
      const double fDelay;
      const art::InputTag fG4ModuleLabel;
      const art::InputTag fOpHitLabel;
      const art::InputTag fOpFlashLabel;
      const double fMinOpHitEnergyFraction;
      mutable std::vector<art::Ptr<sim::OpDetBacktrackerRecord> > priv_OpDetBTRs;
      std::map< art::Ptr < recob::OpFlash >, std::vector < art::Ptr < recob::OpHit > > > fOpFlashToOpHits;


  };//Class
}//namespace

#include "PhotonBackTracker.tcc"

#endif

