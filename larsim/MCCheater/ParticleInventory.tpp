namespace cheat{
  template<typename Evt> //DO NOT USE THIS FUNCTION FROM WITHIN ART! The ParticleInventoryService is designed to impliment these methods as cleanly as possible within the art framework.
    void ParticleInventory::PrepEvent        ( const Evt& evt ){
      fParticleList.clear();
      fMCTruthList.clear();
      fTrackIdToMCTruthIndex.clear();
      this->CheckCanRun(evt);
      if(this->CanRun()){
        this->PrepParticleList(evt);
        this->PrepMCTruthList(evt);
        this->PrepTrackIdToMCTruthIndex(evt);
      }else{throw;}
    }


  template<typename Evt>
    void ParticleInventory::PrepParticleList(const Evt& evt ) const{
      if(!this->CanRun()){throw;}

      if(this->ParticleListReady( )){ //The particle list already exists. Do nothing.
        return;
      }else{ //The particle list needs to be built
        art::Handle<std::vector<simb::MCParticle>> pHandle;
        evt->getByLabel(fG4ModuleLabel, pHandle);

        if(pHandle.failedToGet()){
          /*        mf::LogWarning("ParticleListManager") << "failed to get handle to simb::MCParticle from "
           *                          << fG4ModuleLabel
           *                                            << ", return";
           *                                                              return;*/
          throw;
        }
        for(size_t p=0; p<pHandle->size(); ++p){
          fParticleList.Add( (simb::MCParticle*)(&(pHandle->at(p)))); //Is this still doing a copy? If so, another method should be used.
        }
        fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);
      }
    }

  template<typename Evt>
    void ParticleInventory::PrepTrackIdToMCTruthIndex(Evt evt ) const{
      if(!this->CanRun()){throw;}
      if( this->TrackIdToMCTruthReady( )){//The list already exists. Do nothing.
        return;
      }else{//The List needs to be built
        this->PrepParticleList( evt); //Make sure we have built the particle list for this event
        art::Handle< std::vector<art::Assns<simb::MCParticle,simb::MCTruth, void> > > assnMCParticleTruthHandle;
        evt->getByLabel(fG4ModuleLabel,assnMCParticleTruthHandle);
        if(assnMCParticleTruthHandle.failedToGet()){
          mf::LogWarning("ParticleListManager") << "failed to get handle to art::Assns<simb::MCParticle,simb::MCTruth> from "
            << fG4ModuleLabel
            << ", return";
          return;
        }
        std::cout<<"Size of MCParticleToTruthHandle is: "<<assnMCParticleTruthHandle->size()<<"\n";
        for(unsigned int aN = 0; aN<assnMCParticleTruthHandle->size(); ++aN){
          auto assnsMCParticleToTruth = assnMCParticleTruthHandle->at(aN);
          for( auto assnMCParticleToTruth = assnsMCParticleToTruth.begin(); assnMCParticleToTruth!=assnsMCParticleToTruth.end(); ++assnMCParticleToTruth ){
            art::Ptr<simb::MCParticle> part=assnMCParticleToTruth->first;
            art::Ptr<simb::MCTruth>    mct =assnMCParticleToTruth->second;
            fTrackIdToMCTruthIndex.emplace(part->TrackId(), fMCTruthList.size());
            fMCTruthList.push_back(mct);
          }
        }//end for each Assn<Part, MCT> in event 
      }//End Build TrackIdToTruthList
    }

  template<typename Evt>
    void ParticleInventory::PrepMCTruthList             (Evt& evt ) const{
      if(!this->CanRun()){throw;}
      if(this->MCTruthListReady( ) ){ return;} //If the event is data or if the truth list is already built there is nothing for us to do.
      PrepTrackIdToMCTruthIndex( evt); //TrackIdToMCTruthIndex and MCTruthList are prepared at the same time. The access of information makes this the most convenient way to do so. It is only somewhat more expensive for the memory, but significantly less expensive for time.
    }

  template<typename Evt>
    void ParticleInventory::CheckCanRun( Evt& evt){
      fCanRun= ( !( evt->isRealData() ) ); //If the event is real data, then fCanRun should be false. BackTracking real data makes no sense.
    }
}//end namespace
