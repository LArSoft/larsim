#ifndef TRIGGERALGOMICROBOONE_H
#define TRIGGERALGOMICROBOONE_H

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Event.h"

// Base include
#include "TriggerAlgoBase_service.cc"

// LArSoft include
#include "OpticalDetectorData/OpticalTypes.h"
#include "OpticalDetectorData/PMTTrigger.h"
#include "OpticalDetectorData/FIFOChannel.h"
#include "Simulation/BeamGateInfo.h"
#include "RawData/RawDigit.h"

namespace trigger
{

  /**
     /class TriggerAlgoMicroBoone

     Implementation of abstract TriggerAlgoBase class for MicroBoone experiment.
     
     FillData(art::Event& event) is implemented to read in BeamGateInfo and PMTTrigger array
     from the input data stream.

     Utility function ConvertTime(double time_ns) and ConvertTimeSlice(optdata::TimeSlice_t time_ts) 
     is provided to convert BeamGateInfo and PMTTrigger time units into the time unit used in the
     Trigger logic.

  */
  class TriggerAlgoMicroBoone : public TriggerAlgoBase {
    
  public:
    
    /// Default constructor
    TriggerAlgoMicroBoone(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    /// Default destructor
    virtual ~TriggerAlgoMicroBoone(){};
    
    /// Function to run trigger simulation
    virtual void RunTriggerSim(const art::Event& event){
      FillData(event);
      if(!_sim_done){
	SimTrigger();
	TriggerFIFO(event);
	TriggerWire(event); 
      }
    };

    /// A function to convert input time in nano-seconds into the time unit used in Trigger logic.
    trigdata::TrigTimeSlice_t ConvertTime(double time_ns) const {
      return trigdata::TrigTimeSlice_t( (time_ns*1.e-3) * (_clock_trigger) );
    };

    /**
       A function to convert input time in PMT FEM's timeslice into the time unit used in Trigger logic.
       Note that, in MicroBooNE, we have a frame number and time slice number (though probably it's common to most experiments).
       We convert this into a single time unit.
    */
    trigdata::TrigTimeSlice_t ConvertTimeSlice_PMT(optdata::Frame_t frame, optdata::TimeSlice_t time_ts) const {
      return trigdata::TrigTimeSlice_t( ((double)(time_ts + frame * _pmt_frame_size)) / _clock_pmt * _clock_trigger );
    };

    /**
       A function to convert input time in TPC FEM's timeslice into the time unit used in Trigger logic.
       Note that, in MicroBooNE, we have a frame number and time slice number (though probably it's common to most experiments).
       We convert this into a single time unit.
    */
    trigdata::TrigTimeSlice_t ConvertTimeSlice_TPC(unsigned short frame, unsigned short time_ts) const {
      return trigdata::TrigTimeSlice_t( ((double)(time_ts + frame * _tpc_frame_size)) / _clock_tpc * _clock_trigger );
    };

    /// Alternative to TriggerAlgoBase::IsTriggered(trigdata::TrigTimeSlice_t time) ... take PMTTrigger ptr as input
    bool IsTriggered(const optdata::PMTTrigger* trig) const { 
      return TriggerAlgoBase::IsTriggered( ConvertTimeSlice_PMT(trig->Frame(),trig->TimeSlice()) );
    };

    /// Alternative to TriggerAlgoBase::IsTriggered(trigdata::TrigTimeSlice_t time) ... take FIFOChannel ptr as input
    bool IsTriggered(const optdata::FIFOChannel* fifo) const {
      return TriggerAlgoBase::IsTriggered( ConvertTimeSlice_PMT(fifo->Frame(),fifo->TimeSlice()) );
    };

    /** Alternative to TriggerAlgoBase::IsTriggered(trigdata::TrigTimeSlice_t time) ... take RawDigit ptr as input

       CAUTION: Kazu assumed RawDigit is what is relevant as TPC FEM's raw output data. This data container does not contain
                a frame number as of August 20th 2013. For now, we don't use the frame number, and assume data type is unsigned short.

    */
    bool IsTriggered(const raw::RawDigit* wire) const {
      return TriggerAlgoBase::IsTriggered( ConvertTimeSlice_TPC(0, wire->Samples()) );
    };

    /// Getter function for triggered FIFO channel index array
    const std::set<size_t>* GetIndexes_FIFO() const {return &_readout_index_fifo;};

    /// Getter function for triggered wire index array
    const std::set<size_t>* GetIndexes_Wire() const {return &_readout_index_wire;};

    /**
       Override a function to clear all trigger simulation data.
       The base class clears simulated time window and input trigger timestamps. 
       This class overrides to clear also fifo and wire index std::set.
    */
    virtual void ClearTriggerInfo();

    /// Getter function for TPC frame size
    optdata::TimeSlice_t FrameSizeTPC() const { return _tpc_frame_size;};

    /// Getter function for PMT frame size
    optdata::TimeSlice_t FrameSizePMT() const { return _pmt_frame_size;};

    /// Getter function for Trigger frame size
    optdata::TimeSlice_t FrameSizeTrigger() const { return _trigger_frame_size;};

  protected:

    //
    // Attribute functions
    //

    /// Override implementation of TriggerAlgoBase::FillData()
    virtual void FillData(const art::Event& event);

    /// Variable Initialization ... CAUTION this initializes ALL private attributes (some of which should be set from fhicl file)
    void InitVars();

    /// A function to store a set of index for PMT FIFOChannel that is triggered for readout
    void TriggerFIFO(const art::Event& event);

    /// A function to store a set of index for SimWire that is triggered for readout
    void TriggerWire(const art::Event& event);

    /// A function to read in input fhicl parameters
    void Config(fhicl::ParameterSet const& pset);

    //
    // Attribute variables
    //

    /// A std::set index of FIFO channels that have been triggered for readout
    std::set<size_t> _readout_index_fifo;

    /// A std::set index of wires that have been triggered for readout
    std::set<size_t> _readout_index_wire;

    /// Input data stream module name
    std::string _name_opt_module, _name_gen_module, _name_daq_module;

    /** 
	PMT FEM's sampling frequency in MHz for digitization.
	Used to convert BeamGateInfo time (in ns) into trigdata::TrigTimeSlice_t
    */
    double _clock_pmt;

    /** 
	TPC FEM's sampling frequency in MHz for digitization.
	Used to convert BeamGateInfo time (in ns) into trigdata::TrigTimeSlice_t
    */
    double _clock_tpc;

    /** 
	Trigger module's sampling frequency in MHz for digitization.
	Used to convert BeamGateInfo time (in ns) into trigdata::TrigTimeSlice_t
    */
    double _clock_trigger;

    /// Frame size of PMT FEM readout
    optdata::TimeSlice_t _pmt_frame_size;

    /// Frame size of TPC FEM readout
    optdata::TimeSlice_t _tpc_frame_size;

    /// Frame size of Trigger readout
    optdata::TimeSlice_t _trigger_frame_size;

    /// run, subrun, and event id to keep track of which event has been processed.
    uint32_t _current_run, _current_subrun, _current_event;

    /// Run control booleans
    bool _fifo_sim_done, _wire_sim_done;

  }; /// TriggerAlgoMicroBoone class def

} /// trigger namesapce

DECLARE_ART_SERVICE(trigger::TriggerAlgoMicroBoone, LEGACY)

#endif
