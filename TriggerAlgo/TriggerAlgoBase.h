////////////////////////////////////////////////////////////////////////
// \file TriggerAlgoBase.h
//
// \brief Trigger logic base class. 
//
// \author kazuhiro@nevis.columbia.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef TRIGGERALGOBASE_H
#define TRIGGERALGOBASE_H

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// From this package
#include "TriggerTypes.hh"

// STL
#include <set>

// ROOT
#include <TString.h>

namespace trigger
{
  /**
     Simple read-out trigger logic base class.
     This class applies a logic of data readout trigger.

     Given an input array of readout trigger timestamps, this module applies
     trigger deadtime and readout window to decide which data packet to be 
     readout and stored.

     In the following, I use a notation TS to represent trigdata::TrigTimeSlice_t 
     for short. Note TS is what is used as absolute scale of time in this algorithm.
     Each experiment may implement different reference point or unit to this time
     scale. See TriggerTypes.h for definition.

     std::set<TS> _timestamps should be filled with time-stamps
     (with an arbitrary unit of your choise) of readout candidate triggers.
     This should be done in FillData(art::Event& event) function, which
     is virtual in this header and needs to be implemented in the inherited
     class.

     After filling timestamps, one calls SimTrigger() to analyze and decide
     time windows during which incoming data should be readout and stored.
     This function stores a valid readout time window in the std::map<TS,TS>
     type variable, _time_windows. This map contains the ending time stamp
     of each readout window as a key, and the starting time stamp as a value.

     Once SimTrigger() is called, then one can ask whether a specific time T 
     is within the valid readout windows or not by a function, IsTriggered(TS time).
     This function returns true if the given time stamp is included in any of valid 
     readout window. Else, returns false.

     For further specifications that may have to do with a specific experiment
     or electronics should implement their specific complications in the inherited
     class. This includes FillData() virtual function since different experiments
     certainly have different data format from which readout trigger candidates'
     timestamp is extracted. 

     Note: IsTriggered(TS time) function performs a search of closest readout
     window to "time", and hence could be expensive (though stl::upper_bound is pretty fast).

  */
  class TriggerAlgoBase {

  public:

    /// Default constructor with fhicl parameters
    TriggerAlgoBase(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    /// Default destructor
    virtual ~TriggerAlgoBase(){};

    /// Function to run trigger simulation ... children class may be override
    virtual void RunTriggerSim(const art::Event& event){
      FillData(event);
      SimTrigger();      
    };

    /// Function to clear simulated trigger information
    virtual void ClearTriggerInfo(){
      _sim_done=false;
      _timestamps.clear();
      _time_windows.clear();
    };

    /// Getter for a boolean which "true" value indicates trigger simulation is run already
    bool HasRunTriggerSim() const { return _sim_done;};

    /// Function to check if "time" (input arg.) is within any of valid readout windows or not
    bool IsTriggered(trigdata::TrigTimeSlice_t time) const;

    /// Getter to a const pointer of _time_windows std::map variable
    const std::map<trigdata::TrigTimeSlice_t,trigdata::TrigTimeSlice_t>* GetTimeWindows() const {return &_time_windows;};

    /// Getter to a const pointer of _timestamps std::set variable
    const std::set<trigdata::TrigTimeSlice_t>* GetTriggerTimeStamps() const {return &_timestamps;};


  protected:

    /// Function to fill _timestamps std::set variable ... TO BE IMPLEMENTED in children 
    virtual void FillData(const art::Event& event)=0;

    /// Function to analyze _timestamps and store valid readout windows in _time_windows
    virtual void SimTrigger();

    /// Function to extract fhicl parameters
    void Config(fhicl::ParameterSet const& pset);

    /// stores CANDIDATE readout trigger timestamps
    std::set<trigdata::TrigTimeSlice_t> _timestamps;

    /// stores VALID readout trigger time windows 
    std::map<trigdata::TrigTimeSlice_t,trigdata::TrigTimeSlice_t> _time_windows;

    /// preceeding readout-window from trigger time stamp
    trigdata::TrigTimeSlice_t _preceeding_slices;

    /// proceeding readout-window from trigger time stamp
    trigdata::TrigTimeSlice_t _proceeding_slices;

    /// trigger deadtime AFTER each valid trigger timestamps
    trigdata::TrigTimeSlice_t _deadtime;

    /// run utility boolean, set to true after trigger simulation is run
    bool _sim_done;

  }; // class TriggerAlgoBase
  
} //namespace trigger

#endif
