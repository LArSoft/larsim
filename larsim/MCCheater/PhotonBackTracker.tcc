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

}
