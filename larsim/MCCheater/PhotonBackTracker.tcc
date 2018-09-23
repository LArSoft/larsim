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
      //      if(btrHandle.failedToGet()){
      /*  mf::LogWarning("PhotonBackTracker") << "failed to get handle to     simb::MCParticle from "
       *              << fG4ModuleLabel
       *                          << ", return";*/ //This is now silent as it is expected to    happen every generation run. It is also temporary while we wait for
      /*if( 0 ){ return;} //Insert check for DivRecs here, or don't use validHandle below.
      auto const& divrecHandle = evt.template getValidHandle <std::vector<sim::OpDetDivRec>>(fWavLabel);
      if(divrecHandle.failedToGet()){
        return;
      }*/

      art::fill_ptr_vector(priv_OpDetBTRs, btrHandle);
      //art::fill_ptr_vector(priv_DivRecs, divrecHandle);

      auto compareBTRlambda = [](art::Ptr<sim::OpDetBacktrackerRecord> a, art::Ptr<sim::OpDetBacktrackerRecord> b) {return(a->OpDetNum()<b->OpDetNum());};
      if (!std::is_sorted(priv_OpDetBTRs.begin(),priv_OpDetBTRs.end(),compareBTRlambda)) 
        std::sort(priv_OpDetBTRs.begin(),priv_OpDetBTRs.end(),compareBTRlambda);
      //auto compareDivReclambda = [](art::Ptr<sim::OpDetDivRec> a, art::Ptr<sim::OpDetDivRec> b) {return(a->OpDetNum() < b->OpDetNum());};
      /*if (!std::is_sorted(priv_DivRecs.begin(), priv_DivRecs.end(), compareDivReclambda)) 
        std::sort(priv_DivRecs.begin(), priv_DivRecs.end(), compareDivReclambda);*/
      //art::FindManyP<raw::OpDetWaveform, sim::OpDetDivRec> fp(priv_OpDetBTRs, evt, fWavLabel);// fp;
      //art::FindOneP<raw::OpDetWaveform, sim::OpDetDivRec> fp(priv_OpDetBTRs, evt, fWavLabel);// fp;
      //They come in sorted by BTR. Now make an index matched vector of data_t sorted by BTR. No. I need easy, not efficient. Map of DetNum to data_t. data_t is then channel mapped.
      /*
         if (fp.isValid()){
         for( size_t btr_iter=0; btr_iter<priv_OpDetBTRs.size(); ++btr_iter){
         auto btr=priv_OpDetBTRs.at(btr_iter);
         auto od = btr->OpDetNum();
         auto const& dr = fp.data(btr_iter);
         for(auto& d : dr)
         {
         if(!d) continue;
         priv_od_to_DivRec[od]=*d;//->ref();
         }

         }
         }else{throw cet::exception("PhotonBackTracker")<<"find Waveforms and DivRecs from BTRs failed.";}
         */

      return;
    }

    //----------------------------------------------------------------
    //ToDo: Figure out why I get OpHit* out of here instead of art::Ptr.
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

