
#ifndef TRIGGERALGOMICROBOONE_CC
#define TRIGGERALGOMICROBOONE_CC

#include "TriggerAlgoMicroBoone.h"

namespace trigger{

  TriggerAlgoMicroBoone::TriggerAlgoMicroBoone(fhicl::ParameterSet const& pset,
					       art::ActivityRegistry& reg) 
    : TriggerAlgoBase(pset, reg) 
  {

    InitVars();

    ClearTriggerInfo();

    Config(pset);

  }

  //***************************************************************
  void TriggerAlgoMicroBoone::InitVars() {
  //***************************************************************

    _name_opt_module = _name_gen_module = _name_daq_module = "";

    _clock_pmt = _clock_tpc = _clock_trigger = 0;

    _pmt_frame_size = _tpc_frame_size = 0;

    _current_run = _current_subrun = _current_event = -1;

  }

  //***************************************************************
  void TriggerAlgoMicroBoone::ClearTriggerInfo() {
  //***************************************************************

    TriggerAlgoBase::ClearTriggerInfo();

    _readout_index_fifo.clear();

    _readout_index_wire.clear();

    _current_run = _current_subrun = _current_event = -1;

    _fifo_sim_done = _wire_sim_done = _sim_done = false;

  }  

  //***************************************************************
  void TriggerAlgoMicroBoone::Config(fhicl::ParameterSet const& pset) {
  //***************************************************************

    _name_opt_module     = pset.get< std::string          >("ModuleName_PMT");

    _name_gen_module     = pset.get< std::string          >("ModuleName_BeamSim");

    _name_daq_module     = pset.get< std::string          >("ModuleName_DAQ");

    _pmt_frame_size      = pset.get< optdata::TimeSlice_t >("FrameSize_PMT");

    _tpc_frame_size      = pset.get< optdata::TimeSlice_t >("FrameSize_TPC");

    _trigger_frame_size  = pset.get< optdata::TimeSlice_t >("FrameSize_Trigger");

    _clock_pmt           = pset.get< double               >("ClockFrequency_PMT");

    _clock_tpc           = pset.get< double               >("ClockFrequency_TPC");

    _clock_trigger       = pset.get< double               >("ClockFrequency_Trigger");

  }

  //***************************************************************
  void TriggerAlgoMicroBoone::FillData(const art::Event& event) {
  //***************************************************************

    // Check if data is already filled or not. If filled, return
    if( _current_run == event.run() && _current_event == event.id().event() )
      
      return;

    // Clear variables as this is a new data filling trial
    ClearTriggerInfo();

    //
    // Read in PMTTrigger array ... handle the case if there's no PMTTrigger 
    //
    std::vector<const optdata::PMTTrigger*> pmtTrigArray;
    try {
      event.getView(_name_opt_module, pmtTrigArray);
    }
    catch ( art::Exception const& err ) {
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;
    }

    //
    // Read in BeamGateInfo array ... handle the case if there's no BeamGateInfo 
    //
    std::vector<const sim::BeamGateInfo*> beamGateArray;
    try{
      event.getView(_name_gen_module, beamGateArray);
    }
    catch ( art::Exception const& err ) {
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;
    }

    //
    // Loop over readout trigger candidate time stamps, and fill _timestamps std::set.
    //

    // First, PMTTrigger array
    for(size_t index=0; index < pmtTrigArray.size(); ++index){

      const optdata::PMTTrigger* trig(pmtTrigArray.at(index));

      trigdata::TrigTimeSlice_t timestamp = ConvertTimeSlice_PMT(trig->Frame(), trig->TimeSlice());

      _timestamps.insert(timestamp);

    }

    // Second, BeamGate array
    for(size_t index=0; index < beamGateArray.size(); ++index){

      const sim::BeamGateInfo* trig(beamGateArray.at(index));

      trigdata::TrigTimeSlice_t timestamp = ConvertTime(trig->Start());

      _timestamps.insert(timestamp);

    }

    // Data is filled. Update run, subrun and event id
    _current_run = event.run();
    _current_event = event.id().event();

    // Done.
  }


  //***************************************************************
  void TriggerAlgoMicroBoone::TriggerFIFO(const art::Event& event) {
  //***************************************************************
    
    if(_fifo_sim_done) return;

    _readout_index_fifo.clear();
    
    //
    // Read in FIFOChannel array ... handle the case if there's no FIFOChannel
    //
    std::vector<const optdata::FIFOChannel*> pmtFIFOArray;
    try {
      event.getView(_name_opt_module, pmtFIFOArray);
    }
    catch ( art::Exception const& err ) {
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;
    }
    
    for(size_t index=0; index < pmtFIFOArray.size(); ++index){
      
      const optdata::FIFOChannel* fifo(pmtFIFOArray.at(index));
      
      if(IsTriggered(fifo)) _readout_index_fifo.insert(index);

    }

  }

  //***************************************************************
  void TriggerAlgoMicroBoone::TriggerWire(const art::Event& event) {
  //***************************************************************

    if(_wire_sim_done) return;

    _readout_index_wire.clear();

    //
    // Read in RawDigits array ... handle the ase if there's no RawDigit
    //
    std::vector<const raw::RawDigit*> tpcWireArray;
    try{
      event.getView(_name_daq_module, tpcWireArray);
    }
    catch ( art::Exception const& err ) {
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;
    }

    for(size_t index=0; index < tpcWireArray.size(); ++index){

      const raw::RawDigit* wire(tpcWireArray.at(index));

      if(IsTriggered(wire)) _readout_index_wire.insert(index);

    }

  }
  
} // namespace trigger

DEFINE_ART_SERVICE(trigger::TriggerAlgoMicroBoone)

#endif 
