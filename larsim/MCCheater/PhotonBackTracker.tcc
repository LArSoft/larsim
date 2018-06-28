namespace cheat{


  //----------------------------------------------------------------
  template<typename Evt>
    const bool PhotonBackTracker::CanRun(Evt const& evt) {
      return ! ( evt.isRealData() ) ;
    }


  //----------------------------------------------------------------
  template<typename Evt>
    void PhotonBackTracker::PrepOpDetBTRs(Evt const& evt)
    {
      if(this->BTRsReady()){ return;}
      auto const& btrHandle = evt.template getValidHandle < std::vector < sim::OpDetBacktrackerRecord > > (fG4ModuleLabel);
      if(btrHandle.failedToGet()){
        /*  mf::LogWarning("PhotonBackTracker") << "failed to get handle to     simb::MCParticle from "
         *              << fG4ModuleLabel
         *                          << ", return";*/ //This is now silent as it is expected to    happen every generation run. It is also temporary while we wait for
        return;
      }

      art::fill_ptr_vector(priv_OpDetBTRs, btrHandle);

      auto compareBTRlambda = [](art::Ptr<sim::OpDetBacktrackerRecord> a, art::Ptr<sim::OpDetBacktrackerRecord> b) {return(a->OpDetNum()<b->OpDetNum());};
      if (!std::is_sorted(priv_OpDetBTRs.begin(),priv_OpDetBTRs.end(),compareBTRlambda)) std::sort(priv_OpDetBTRs.begin(),priv_OpDetBTRs.end(),compareBTRlambda);

      return;
    }

  //----------------------------------------------------------------
  template<typename Evt>
    void PhotonBackTracker::PrepOpFlashToOpHits( Evt const& evt)
    {
      if(this->OpFlashToOpHitsReady()){ return;}
      auto const& flashHandle = evt.template getValidHandle < std::vector < recob::OpFlash > > (fOpFlashLabel);
      std::vector<art::Ptr<recob::OpFlash>> tmpVec;
      art::fill_ptr_vector(tmpVec, flashHandle);

      if(flashHandle.failedToGet()){
        mf::LogWarning("PhotonBackTracker")<<" failed to get handle to recob::OpFlash. Has reco run yet?";
        return;
      }
      auto const& fmp = art::FindManyP<recob::OpHit>(flashHandle, evt, fOpHitLabel.label());
      for ( size_t i=0; i<flashHandle->size(); i++)
        fOpFlashToOpHits[tmpVec.at(i)] = fmp.at(i);
    }

  //----------------------------------------------------------------
  template<typename Evt>
    void PhotonBackTracker::PrepEvent( Evt const& evt)
    {
      if( !(this->CanRun( evt ) ) ){
        throw cet::exception("PhotonBackTracker")
          <<"PhotonBackTracker cannot function."
          <<"Is this file real data?";
      }
      priv_OpDetBTRs.clear();
      this->PrepOpDetBTRs(evt);
      this->PrepOpFlashToOpHits(evt);
    } 



}

