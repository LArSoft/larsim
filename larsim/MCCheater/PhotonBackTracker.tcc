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
      FindOneP<raw::OpDetWaveform, sim::OpDetDivRec>(priv_OpDetBTRs, evt, fWavLabel) fp;
      //They come in sorted by BTR. Now make an index matched vector of data_t sorted by BTR. No. I need easy, not efficient. Map of DetNum to data_t. data_t is then channel mapped.
      if (fp.isValid()){
        for( auto& btr : priv_OpDetBTRs){
          //BTRs are now sorted, but I want a map anyway. I can do better, but this is quick and robust.
          auto check = od_to_chanDiv.emplace(std::make_pair(btr.OpDetNum(), fp.data(btr)));
          if(! check.second){
            throw cet::exception("PhotonBackTracker")<<"Trying to add a ChannelDivRec by BTR for an already existing entry.";
          }

        }
      }else{throw cet::exception("PhotonBackTracker")<<"find Waveforms and DivRecs from BTRs failed.";}
      
      return;
    }

  //----------------------------------------------------------------
  template<typename Evt>
    void PhotonBackTracker::PrepOpFlashToOpHits( Evt const& evt)
    {
      if(this->OpFlashToOpHitsReady()){ return;}
      auto const& flashHandle = evt.template getValidHandle < std::vector < recob::OpFlash > > (fOpFlashLabel.label());
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

