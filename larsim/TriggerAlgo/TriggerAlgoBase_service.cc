////////////////////////////////////////////////////////////////////////
//
//  \file TriggerAlgoBase_service.cc
//
////////////////////////////////////////////////////////////////////////

#ifndef TRIGGERALGOBASE_CC
#define TRIGGERALGOBASE_CC

#include "TriggerAlgoBase.h"

namespace trigger{

  //****************************************************************************
  TriggerAlgoBase::TriggerAlgoBase(fhicl::ParameterSet const& pset, art::ActivityRegistry& /* reg */) {
  //****************************************************************************

    ClearTriggerInfo();

    Config(pset);

  }

  //****************************************************************************
  void TriggerAlgoBase::Config(fhicl::ParameterSet const& pset) {
  //****************************************************************************
    
    _preceeding_slices = pset.get< int >("PreceedingWindow");

    _proceeding_slices = pset.get< int >("ProceedingWindow");

    _deadtime          = pset.get< int >("DeadTime");

  }


  //****************************************************************************
  void TriggerAlgoBase::SimTrigger() {
  //****************************************************************************

    if(_sim_done) return;

    _time_windows.clear();

    trigdata::TrigTimeSlice_t last_timestamp=0;

    trigdata::TrigTimeSlice_t window_begin=0;

    trigdata::TrigTimeSlice_t window_end=0;

    for(std::set<trigdata::TrigTimeSlice_t>::const_iterator iter(_timestamps.begin());
	iter != _timestamps.end();
	++iter){
      
      if(!(last_timestamp) || (*iter) > (last_timestamp + _deadtime)) {
	
	window_begin = ((*iter) > _preceeding_slices) ? ((*iter) - _preceeding_slices) : 0;
	
	window_end   = (*iter) + _proceeding_slices;

	_time_windows.insert(std::make_pair(window_end,window_begin));
	
	last_timestamp=(*iter);
	
      }

    }

    _sim_done = true;

  }

  //****************************************************************************
  bool TriggerAlgoBase::IsTriggered(trigdata::TrigTimeSlice_t time) const {
  //****************************************************************************
    
    if(!_time_windows.size()) return false;

    std::map<trigdata::TrigTimeSlice_t,trigdata::TrigTimeSlice_t>::const_iterator start_time(_time_windows.lower_bound(time));

    if(start_time==_time_windows.end()) return false;

    else return ((*start_time).second < time);
    
  }
  
  
  //DEFINE_ART_SERVICE(TriggerAlgoBase)

} // namespace trigger

#endif 
