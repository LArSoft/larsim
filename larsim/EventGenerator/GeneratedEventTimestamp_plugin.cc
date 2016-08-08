/**
 * @file    GeneratedEventTimestamp_plugin.cc
 * @brief   Assigns an empty event a time stamp from the clock
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * 
 * This file defines a plug in for art framework; no interface is needed since
 * users will interface with the abstract base class.
 */

// C/C++ standard libraries
#include <chrono>
#include <random>

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "art/Framework/Core/EmptyEventTimestampPlugin.h"

// Event generation namespace
namespace evgen {
  
  /**
   * @brief Plugin to assign an empty event a time stamp from the clock
   * @see art::EmptyEventTimestampPlugin
   *
   * The plug in returns a time stamp that is taken from the current time 
   * on the execution node, in nanoseconds.
   */
  class GeneratedEventTimestamp: public art::EmptyEventTimestampPlugin {
      public:
    
    /// Constructor: nothing specific
    GeneratedEventTimestamp(fhicl::ParameterSet const& pset);
    
    
    /// Returns the time stamp for the specified event
    virtual art::Timestamp eventTimestamp(art::EventID const&) override;
    
    /// Resets the status; since this plug in is stateless, this is a no-op
    virtual void rewind() override {}
    
    
  }; // class GeneratedEventTimestamp
  
  
  
  
} // namespace evgen

//------------------------------------------------------------------------------
//--- Implementation
//---
//---

namespace evgen {
  namespace details {
    
    /// Class reading a clock and converting the value to a specific unit
    template <typename Clock, typename Duration>
    class TimeInUnitsBase {
        public:
      
      long long int operator() () { return read_clock(); }
      
      static long long int read_clock()
        {
          // obtain from the high resolution clock the current time,
          // from the "epoch":
          typename Clock::duration now = Clock::now().time_since_epoch();
          
          // convert the time in the requested unit
          return std::chrono::duration_cast<Duration>(now).count();
        } // operator()
      
    }; // class TimeInUnitsBase


    /// Class reading a clock and converting the value to a specific unit;
    /// if the unit is more precise than the clock, random padding fills the gap
    template <typename Clock, typename Duration, typename = void>
    class TimeInUnits: public TimeInUnitsBase<Clock, Duration> {
        public:
      using Base_t = TimeInUnitsBase<Clock, Duration>;
      
      long long int operator() () { return Base_t::read_clock(); }
    }; // class TimeInUnits<>


    // Implementation of the random-gap-filling version
    template <typename Clock, typename Duration>
    class TimeInUnits<Clock, Duration, typename std::enable_if<
        (Clock::period::num * Duration::period::den > Duration::period::num * Clock::period::den)
        >::type
      >
      : public TimeInUnitsBase<Clock, Duration>
    {
      using Base_t = TimeInUnitsBase<Clock, Duration>;
      
      // if the period of the clock is larger than the unit we are requested,
      // there will be some padding
      using ClockPeriod = typename Clock::period;
      using ReqPeriod = typename Duration::period;
      // requested clock / unit:
      using PeriodRatio = std::ratio<
        ClockPeriod::num * ReqPeriod::den, ReqPeriod::num * ClockPeriod::den
        >;
      static constexpr intmax_t paddingfactor
        = PeriodRatio::num / PeriodRatio::den; // > 1 enforced by enable_if<>
      
        public:
      TimeInUnits(): engine(), flat(0, paddingfactor - 1) {}
      
      /// Return the clock value with random padding added
      long long int operator() ()
        { return Base_t::read_clock() + flat(engine); }
      
        private:
      std::default_random_engine engine;
      std::uniform_int_distribution<long long int> flat;
    }; // class TimeInUnits<> (padded)
  
  } // namespace details
} // namespace evgen


//------------------------------------------------------------------------------
evgen::GeneratedEventTimestamp::GeneratedEventTimestamp
  (fhicl::ParameterSet const& pset)
  : art::EmptyEventTimestampPlugin(pset)
{
  mf::LogInfo("GeneratedEventTimestamp")
    << "Timestamp plugin: timestamp from local clock time in nanoseconds.";
} // evgen::GeneratedEventTimestamp::GeneratedEventTimestamp()


//------------------------------------------------------------------------------
art::Timestamp evgen::GeneratedEventTimestamp::eventTimestamp
  (art::EventID const&)
{
  // obtain from the high resolution clock the current time, from the "epoch",
  // in nanoseconds; if the clock is less precise than the nanosecond,
  // the precision gap is filled with randomness
  details::TimeInUnitsBase
    <std::chrono::high_resolution_clock, std::chrono::nanoseconds>
    get_time;
  
  const long long int now_ns = get_time();
  
  // convert into a timestamp
  art::Timestamp ts(now_ns);
  
  LOG_DEBUG("GeneratedEventTimestamp") << "Generated time stamp: "
    << ts.value();
  
  return ts;
} // evgen::GeneratedEventTimestamp::eventTimestamp()

// make art aware that we have a plugin
DEFINE_ART_EMPTYEVENTTIMESTAMP_PLUGIN(evgen::GeneratedEventTimestamp)

