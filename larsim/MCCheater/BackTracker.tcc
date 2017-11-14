namespace cheat{

  //--------------------------------------------------------------------
  template<typename Evt> //DO NOT USE THIS FUNCTION FROM WITHIN ART! The BackTrackerService is designed to impliment these methods as cleanly as possible within the art framework. This is intended for gallery users.
    void BackTracker::PrepEvent        (const Evt& evt ){
      if( !( this->CanRun( evt ) ) ){       
        throw cet::exception("ParticleInventory")
          << "Particle Inventory cannot function. "
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
      art::fill_ptr_vector(fSimChannels, simChannelsHandle); 
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
