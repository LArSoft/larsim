////////////////////////////////////////////////////////////////////////////
//
// \file BackTracker.h
// \brief Functions needed by the BackTracker service in order to connect truth information with reconstruction.
//
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////
#ifndef CHEAT_BACKTRACKER_H
#define CHEAT_BACKTRACKER_H

#include <vector>

#include "fhiclcpp/types/Atom.h"

#include "larsim/MCCheater/ParticleInventory.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace fhicl { class ParameterSet; }
namespace detinfo { class DetectorClocks; }
namespace geo { class GeometryCore; }
namespace recob { class SpacePoint; }

namespace cheat{

  class BackTracker{
    public:

      //Structure for configuration parameters. (fhicl validation)
      struct fhiclConfig{
        fhicl::Atom<art::InputTag> G4ModuleLabel{fhicl::Name("G4ModuleLabel"), fhicl::Comment("The label of the LArG4   module used to produce the art file we will be using."), "largeant"};
        fhicl::Atom<art::InputTag> SimChannelModuleLabel{fhicl::Name("SimChannelModuleLabel"), fhicl::Comment("The label of the module containing the sim::SimChannel product."), G4ModuleLabel()}; // -- D.R. label not required, if not provided defaults to the value of G4ModuleLabel
        fhicl::Atom<art::InputTag> DefaultHitModuleLabel{fhicl::Name("DefaultHitModuleLabel"), fhicl::Comment("The label  of the module used to produce the hits in the art file we will default to when no hitlist is provided."), "hitfd"};
        fhicl::Atom<double> MinHitEnergyFraction{fhicl::Name("MinHitEnergyFraction"), fhicl::Comment("The minimum     contribution an energy deposit must make to a Hit to be considered part of that hit."),0.010};
        fhicl::Atom<bool> OverrideRealData{fhicl::Name("OverrideRealData"),fhicl::Comment("Option when overlaying simulation on real data, to tell the backtracker to continue even if event looks like data."),false};
      };


      BackTracker(const fhiclConfig& config,
          const cheat::ParticleInventory* partInv,
          const geo::GeometryCore* geom,
          const detinfo::DetectorClocks* detClock );
      BackTracker(const fhicl::ParameterSet& pSet,
          const cheat::ParticleInventory* partInv,
          const geo::GeometryCore* geom,
          const detinfo::DetectorClocks* detClock );
      //I may need to include this to delete copy of service providers.
      BackTracker(BackTracker const&) = delete;

      template<typename Evt>
        void PrepEvent ( const Evt& evt );

      template<typename Evt>
        void PrepSimChannels ( const Evt& evt );

//      template<typename Evt>
//        void PrepAllHitList ( const Evt& evt);

      //-----------------------------------------------------
      template<typename Evt>
        bool CanRun(const Evt& evt){ return ( !(evt.isRealData()) || fOverrideRealData);}

      //-----------------------------------------------------
      template<typename Evt>
        const std::vector< art::Ptr< recob::Hit > > SpacePointToHits_Ps(art::Ptr< recob::SpacePoint> const& spt, const Evt& evt) const;

      //-----------------------------------------------------
      template<typename Evt>
        const std::vector<  double > SpacePointToXYZ( art::Ptr< recob::SpacePoint > const& spt, const Evt& evt) const;


      //-----------------------------------------------------
      void ClearEvent();


      bool SimChannelsReady() const { return !( fSimChannels.empty() ); }
//      bool AllHitListReady () const { return !( fAllHitList.empty() ); }

      const std::vector<art::Ptr<sim::SimChannel>>& SimChannels() const { return fSimChannels; }
      //All Hit List would go here. We explicitly choose not to include it, as the user should not be using backtracker to access Hits. This could change in a concievable future use case where we also allow the user to define what the "AllHitList" should be, though this would have ramifications on other functions.

      std::vector<const sim::IDE* >   TrackIdToSimIDEs_Ps(int const& id) const;
      std::vector<const sim::IDE* >   TrackIdToSimIDEs_Ps(int const& id, const geo::View_t view) const;

      art::Ptr<sim::SimChannel> FindSimChannel( raw::ChannelID_t channel ) const;

      std::vector< sim::TrackIDE > ChannelToTrackIDEs(raw::ChannelID_t channel, const double hit_start_time, const double hit_end_time) const;


      //Track IDEs cannot be returned as pointers, as they dont exist in the data product, and we will not be storing them.
      std::vector< sim::TrackIDE> HitToTrackIDEs(recob::Hit const& hit) const;
      std::vector< sim::TrackIDE> HitToTrackIDEs(art::Ptr<recob::Hit> const& hit) const { return this->HitToTrackIDEs(*hit);}

      std::vector< int > HitToTrackIds(recob::Hit const& hit) const ;
      //   std::vector< const int> HitToTrackId(art::Ptr<recob::Hit> const& hit) { return this->HitToTrackId(*hit); }

      std::vector<sim::TrackIDE> HitToEveTrackIDEs(recob::Hit const& hit) const;
      std::vector<sim::TrackIDE> HitToEveTrackIDEs(art::Ptr<recob::Hit> const& hit) const{ return this->HitToEveTrackIDEs(*hit);}

      //I will not return these by copy, as that could get very large very quickly.
      std::vector< art::Ptr<recob::Hit> > TrackIdToHits_Ps( const int& tkId, std::vector< art::Ptr< recob::Hit > > const& hitsIn ) const;
//I am not allowing this function for now, as caching the allhitlist indiscriminately is a catastrophically bad idea.
//      std::vector< art::Ptr<recob::Hit> > TrackIdToHits_Ps( const int& tkId ) const
//      {return this->TrackIdToHits_Ps(tkId, fAllHitList); }

      std::vector< std::vector< art::Ptr<recob::Hit> > > TrackIdsToHits_Ps( std::vector<int> const& tkIds, std::vector< art::Ptr< recob::Hit > > const& hitsIn ) const;
//      std::vector< std::vector< art::Ptr<recob::Hit> > > TrackIdsToHits_Ps( std::vector<int> const& tkIds ) const
//      {return this->TrackIdsToHits_Ps(tkIds, fAllHitList);}

      std::vector< sim::IDE > HitToAvgSimIDEs ( recob::Hit const& hit) const;
      std::vector< sim::IDE > HitToAvgSimIDEs ( art::Ptr<recob::Hit> hit) const{ return this->HitToAvgSimIDEs(*hit);}

      std::vector< const sim::IDE* > HitToSimIDEs_Ps (recob::Hit const& hit) const;
      std::vector< const sim::IDE* > HitToSimIDEs_Ps (art::Ptr< recob::Hit > const& hit) const { return this->HitToSimIDEs_Ps (*hit); }

      //   std::vector< sim::IDE > HitToSimIDEs (recob::Hit const& hit);
      //   std::vector< const sim::IDE > HitToSimIDEs (art::Ptr< recob::Hit > const& hit) { return this->HitToSimIDEsPs (*hit); }

      std::vector<double> SimIDEsToXYZ( std::vector< sim::IDE > const& ides) const;
      std::vector<double> SimIDEsToXYZ( std::vector< const sim::IDE* > const& ide_Ps) const;

      std::vector<double> HitToXYZ(const recob::Hit& hit) const;
      std::vector<double> HitToXYZ(art::Ptr<recob::Hit> const& hit) const{ return this->HitToXYZ(*hit);}



      double HitCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const& hits) const;
      double HitChargeCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const& hits) const;

      double HitCollectionEfficiency( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const& hits, std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view) const;


      double HitChargeCollectionEfficiency( std::set<int> trackIds, std::vector< art::Ptr<recob::Hit> > const& hits,        std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view) const;

      std::set<int> GetSetOfTrackIds() const { return fPartInv->GetSetOfTrackIds();}
      std::set<int> GetSetOfEveIds() const { return fPartInv->GetSetOfEveIds();}

      std::set<int> GetSetOfTrackIds( std::vector< art::Ptr< recob::Hit > > const& hits ) const;
      std::set<int> GetSetOfEveIds( std::vector< art::Ptr< recob::Hit > > const& hits ) const;

      std::vector<  double> SpacePointHitsToWeightedXYZ(std::vector<art::Ptr<recob::Hit>> const& hits) const;

    private:
      const cheat::ParticleInventory* fPartInv; //The constructor needs to put something in here
      const geo::GeometryCore* fGeom;
      const detinfo::DetectorClocks* fDetClocks;
      const art::InputTag       fG4ModuleLabel;
      const art::InputTag       fSimChannelModuleLabel;
      const art::InputTag       fHitLabel;
      const double              fMinHitEnergyFraction;
      const bool                fOverrideRealData;

      mutable std::vector<art::Ptr<sim::SimChannel>>       fSimChannels;
//      mutable std::vector< art::Ptr<recob::Hit> >          fAllHitList;

  };//end class BackTracker

}//end namespace cheat

#include "BackTracker.tcc"

#endif //CHEAT_BACKTRACKER_H
