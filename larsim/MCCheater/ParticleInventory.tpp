namespace cheat{

  //--------------------------------------------------------------------
  template<typename Evt> //DO NOT USE THIS FUNCTION FROM WITHIN ART! The ParticleInventoryService is designed to impliment these methods as cleanly as possible within the art framework.
    void ParticleInventory::PrepEvent        (const Evt& evt ){
      if(this->CanRun(evt)){ 
        throw cet::exception("ParticleInventory") 
        << "Particle Inventory cannot function. "
        << "Is this file real data?";
      }
      fParticleList.clear();
      fMCTObj.fMCTruthList.clear();
      fMCTObj.fTrackIdToMCTruthIndex.clear();
      this->PrepParticleList(evt);
      this->PrepMCTruthList(evt);
      this->PrepTrackIdToMCTruthIndex(evt);
    }


  //--------------------------------------------------------------------
  template<typename Evt>
    void ParticleInventory::PrepParticleList(const Evt& evt ) const{

      if(this->ParticleListReady( )){ //The particle list already exists. Do nothing.
        return;
      }
      //The particle list needs to be built
      //We use auto so that we(the compiler) can determine which type we need for either art or gallery.
      const auto& partVecIn = *(evt.getValidHandle< std::vector< simb::MCParticle > > (fG4ModuleLabel));


      for(const auto& partIn : partVecIn){
        fParticleList.Add(new simb::MCParticle(partIn)); //Is this still doing a copy? If so, another method should be used.
      }
      fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);
    }

  //--------------------------------------------------------------------
  template<typename Evt> //I may want to make this function private.
    void ParticleInventory::PrepMCTruthListAndTrackIdToMCTruthIndex(const Evt& evt ) const{
      if( this->TrackIdToMCTruthReady() && this->MCTruthListReady( ) ){ return;} 
      this->PrepParticleList( evt); //Make sure we have built the particle list for this event
      const auto& mcpmctAssnVecIn = *( evt.getValidHandle< art::Assns<simb::MCParticle,simb::MCTruth> >(fG4ModuleLabel));
      //std::cout<<"Size of MCParticleToTruthHandle is: "<<assnMCParticleTruthHandle->size()<<"\n";
      for( const auto& mcpmctAssnIn : mcpmctAssnVecIn ){
        const art::Ptr<simb::MCParticle>& part=mcpmctAssnIn.first;
        const art::Ptr<simb::MCTruth>&    mct =mcpmctAssnIn.second;
        fMCTObj.fTrackIdToMCTruthIndex.emplace(part->TrackId(), fMCTObj.fMCTruthList.size());
        fMCTObj.fMCTruthList.push_back(mct);
      }
    }

  template<typename Evt>
    void ParticleInventory::PrepMCTruthList             (const Evt& evt ) const{
      if(this->MCTruthListReady( ) ){ return;} //If the event is data or if the truth list is already built there is nothing for us to do.
      PrepMCTruthListAndTrackIdToMCTruthIndex( evt); //TrackIdToMCTruthIndex and MCTruthList are prepared at the same time. The access of information makes this the most convenient way to do so. It is only somewhat more expensive for the memory, but significantly less expensive for time.
    }

  template<typename Evt>
    void ParticleInventory::PrepTrackIdToMCTruthIndex(const Evt& evt ) const{
      if(this->TrackIdToMCTruthReady()){ return;} //The list already exists. Do nothing.
      PrepMCTruthListAndTrackIdToMCTruthIndex( evt); //TrackIdToMCTruthIndex and MCTruthList are prepared at the same time. The access of information makes this the most convenient way to do so. It is only somewhat more expensive for the memory, but significantly less expensive for time.
    }

  template<typename Evt>
    bool ParticleInventory::CanRun(const Evt& evt) const{
      return !(evt.isRealData());
    }

}//end namespace
