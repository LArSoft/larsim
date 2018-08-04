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
      //auto const& flashHandle = evt.template getValidHandle < std::vector < recob::OpFlash > > (fOpFlashLabel);
      art::Handle< std::vector<recob::OpFlash> > flashHandle;
      //= evt.template getValidHandle < std::vector < recob::OpFlash > > (fOpFlashLabel);
      evt.template getByLabel(fOpFlashLabel, flashHandle);
      std::vector<art::Ptr<recob::OpFlash>> tmpVec;
      art::fill_ptr_vector(tmpVec, flashHandle);
      art::Handle< art::Assns<recob::OpFlash,recob::OpHit> > assnFlashToOpHit;
      evt.getByLabel(fOpFlashLabel,assnFlashToOpHit);

      //if(assnFlashToOpHit.failedToGet()){
/*      if(flashHandle.failedToGet()){
        mf::LogWarning("PhotonBackTracker")<<" failed to get handle to recob::OpFlash. Has reco run yet?";
        return;
      }*/

//      std::vector< const std::< art::Ptr<recob::OpHit> > 
      auto tmp = util::GetAssociatedVectorManyP(assnFlashToOpHit, flashHandle);

      for ( size_t i=0; i<flashHandle->size(); i++){
        //auto shape1= tmpVec.at(i);
        auto check = fOpFlashToOpHits.emplace(tmpVec.at(i),tmp.at(i)); //.[tmpVec.at(i)] = fmp.at(i);
        if(check.second==false){mf::LogWarning("PhotonBackTracker")<<"Failed to insert Flash,vec<OpHit> record. Do they exist?";}
      }
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

