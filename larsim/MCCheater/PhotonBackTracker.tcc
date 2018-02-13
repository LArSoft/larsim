namespace cheat{


  //----------------------------------------------------------------
  template<typename Evt>
    const bool PhotonBackTracker::CanRun(const Evt& evt) {
      return ! ( evt.isRealData() ) ;
    }


  //----------------------------------------------------------------
  template<typename Evt>
    void PhotonBackTracker::PrepOpDetBTRs(const Evt& evt)
    {
      if(this->BTRsReady()){ return;}
      const auto& btrHandle = evt.template getValidHandle < std::vector < sim::OpDetBacktrackerRecord > > (fG4ModuleLabel);
      if(btrHandle.failedToGet()){
        /*  mf::LogWarning("PhotonBackTracker") << "failed to get handle to     simb::MCParticle from "
         *              << fG4ModuleLabel
         *                          << ", return";*/ //This is now silent as it is expected to    happen every generation run. It is also temporary while we wait for
      }
    }

  //----------------------------------------------------------------
  template<typename Evt>
  void PhotonBackTracker::PrepEvent( const Evt& evt)
  {
    if( !(this->CanRun( evt ) ) ){
      throw cet::exception("PhotonBackTracker")
        <<"PhotonBackTracker cannot function."
        <<"Is this file real data?";
    }
    priv_OpDetBTRs.clear();
    this->PrepOpDetBTRs(evt);
  } 

  //----------------------------------------------------- /*NEW*/
  template<typename Evt>
    const std::vector<art::Ptr<recob::OpHit>> PhotonBackTracker::OpFlashToOpHits_Ps(art::Ptr<recob::OpFlash>& flash_P, const Evt& evt) const
    {//There is not "non-pointer" version of this because the art::Ptr is needed to look up the assn
      //std::vector<art::Ptr< recob::OpHit > > hits_Ps = FindManyP(flash_P, evt, opHitLabel);
      return FindManyP(flash_P, evt, fOpHitLabel);

    }

  //----------------------------------------------------- /*NEW*/
  template<typename Evt>
    const std::vector<double> PhotonBackTracker::OpFlashToXYZ(art::Ptr<recob::OpFlash>& flash_P, Evt& evt) const
    {
      return this->OpHitsToXYZ( this->OpFlashToOpHits_Ps(flash_P, evt)  );
    }

  //----------------------------------------------------- /*NEW*/
  template<typename Evt>
    const std::set<int> PhotonBackTracker::OpFlashToTrackIds(art::Ptr<recob::OpFlash>& flash_P, Evt& evt) const{
      std::vector<art::Ptr<recob::OpHit> > opHits_Ps = this->OpFlashToOpHits_Ps(flash_P, evt);
      std::set<int> ids;
      for( auto& opHit_P : opHits_Ps){
        for( int& id : this->OpHitToTrackIds(opHit_P) ){
          ids.insert( id) ;
        } // end for ids
      }// end for opHits
      return ids;
    }// end OpFlashToTrackIds

}

