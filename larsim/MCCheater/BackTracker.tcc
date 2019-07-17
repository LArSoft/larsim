////////////////////////////////////////////////////////////////////////////
//
// \file BackTracker.tcc
// \brief Template functions needed by the BackTracker service in order to connect truth information with reconstruction. Templates are for functions that require the event (as the "event" is different for art and Gallery).
//
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/SpacePoint.h"

namespace cheat{

  //--------------------------------------------------------------------
  template<typename Evt> //DO NOT USE THIS FUNCTION FROM WITHIN ART! The BackTrackerService is designed to impliment these methods as cleanly as possible within the art framework. This is intended for gallery users.
    void BackTracker::PrepEvent        (const Evt& evt ){
      if( !( this->CanRun( evt ) ) ){       
        throw cet::exception("BackTracker")
          << "BackTracker cannot function. "
          << "Is this file real data?";
      }
      fSimChannels.clear();
      //      fAllHitList.clear();
      this->PrepSimChannels( evt );
      //this->PrepAllHitList ( evt ); //This line temporarily commented out until I figure out how I want PrepAllHitList to work.

    }

  //--------------------------------------------------------------------
  template<typename Evt>
    void BackTracker::PrepSimChannels (const Evt& evt){
      if(this->SimChannelsReady()){ return;}
      //The SimChannels list needs to be built.
      const auto& simChannelsHandle = evt.template getValidHandle<std::vector<sim::SimChannel>>(fG4ModuleLabel);
      //failedToGet for a valid handle will always be false.
//      if(simChannelsHandle.failedToGet()){
//        /*  mf::LogWarning("BackTracker") << "failed to get handle to simb::MCParticle from "
//            << fG4ModuleLabel
//            << ", return";*/ //This is now silent as it is expected to happen every generation run. It is also temporary while we wait for
//        return;
//      }

      art::fill_ptr_vector(fSimChannels, simChannelsHandle); 

      auto comparesclambda = [](art::Ptr<sim::SimChannel> a, art::Ptr<sim::SimChannel> b) {return(a->Channel()<b->Channel());};
      if (!std::is_sorted(fSimChannels.begin(),fSimChannels.end(),comparesclambda)) std::sort(fSimChannels.begin(),fSimChannels.end(),comparesclambda);

      return;

    }

  //--------------------------------------------------------------------
  /*  template<typename Evt>
      void BackTracker::PrepAllHitList( const Evt& evt){
      if(this->AllHitListReady()){return;}
      const auto& allHitsHandle = evt.template getValidHandle<std::vector<recob::Hit>>(fHitLabel);
      art::fill_ptr_vector(fAllHitList, allHitsHandle);
      }
      */
  //--------------------------------------------------------------------
  template<typename Evt>
    const std::vector< art::Ptr< recob::Hit > > BackTracker::SpacePointToHits_Ps(art::Ptr<recob::SpacePoint> const& spt, const Evt& evt) const{
      std::vector<art::Ptr<recob::SpacePoint>> spv; //This method needs to be rethought. For now I am directly implimenting it as found in the previous backtracker.
      spv.push_back(spt);
      art::FindManyP<recob::Hit> fmh(spv, evt, fHitLabel);
      std::vector< art::Ptr<recob::Hit> > hitv = fmh.at(0);
      return hitv;
    }

  //--------------------------------------------------------------------
  template<typename Evt>
    const std::vector< double > BackTracker::SpacePointToXYZ( art::Ptr< recob::SpacePoint > const& spt, const Evt& evt ) const{
      std::vector<art::Ptr<recob::Hit>> hits=this->SpacePointToHits_Ps(spt, evt);
      return this->SpacePointHitsToWeightedXYZ(hits);
    }

}//end namespace
