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
    void PhotonBackTracker::PrepEvent( Evt const& evt)
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
    const std::vector<art::Ptr<recob::OpHit>> PhotonBackTracker::OpFlashToOpHits_Ps(art::Ptr<recob::OpFlash>& flash_P, Evt const& evt) const
    {//There is not "non-pointer" version of this because the art::Ptr is needed to look up the assn. One could loop the Ptrs and dereference them, but I will not encourage the behavior by building the tool to do it.
      art::FindManyP< recob::OpHit > fmoh(std::vector<art::Ptr<recob::OpFlash>>({flash_P}), evt, fOpHitLabel.label());
      std::vector<art::Ptr<recob::OpHit>> const& hits_Ps = fmoh.at(0);
      //std::vector<art::Ptr< recob::OpHit > > hits_Ps = art::FindManyP< recob::OpHit > fmoh(flash_P, evt, fOpHitLabel);
      return hits_Ps;// art::FindManyP< recob::OpHit >(flash_P, evt, fOpHitLabel);

    }

  //----------------------------------------------------- /*NEW*/
  template<typename Evt>
    const std::vector<double> PhotonBackTracker::OpFlashToXYZ(art::Ptr<recob::OpFlash>& flash_P, Evt& evt) const
    {
      const std::vector< art::Ptr<recob::OpHit>> opHits_Ps = this->OpFlashToOpHits_Ps(flash_P, evt);
      const std::vector<double> retVec = this->OpHitsToXYZ(opHits_Ps);
      return retVec;
    }

  //----------------------------------------------------- /*NEW*/
  template<typename Evt>
    const std::set<int> PhotonBackTracker::OpFlashToTrackIds(art::Ptr<recob::OpFlash>& flash_P, Evt& evt) const{
      std::vector<art::Ptr<recob::OpHit> > opHits_Ps = this->OpFlashToOpHits_Ps(flash_P, evt);
      std::set<int> ids;
      for( auto& opHit_P : opHits_Ps){
        for( const int& id : this->OpHitToTrackIds(opHit_P) ){
          ids.insert( id) ;
        } // end for ids
      }// end for opHits
      return ids;
    }// end OpFlashToTrackIds

}

