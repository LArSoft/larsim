/**
 * @file   EventSeedInputData.h
 * @brief  A data object holding enough data to define a event seed
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   February 18th, 2015
 */


#ifndef LARSIM_RANDOMUTILS_PROVIDERS_EVENTSEEDINPUTDATA_H
#define LARSIM_RANDOMUTILS_PROVIDERS_EVENTSEEDINPUTDATA_H 1

// C/C++ standard libraries
#include <cstdint> // std::uint32_t
#include <string>


namespace sim {
  namespace SeedServiceHelper {
    
    /// Simple data structure with data needed to extract a seed from a event
    class EventSeedInputData {
        public:
      using RunNumber_t    = std::uint32_t;
      using SubRunNumber_t = std::uint32_t;
      using EventNumber_t  = std::uint32_t;
      using TimeValue_t    = std::uint64_t;
      
      /// @{
      /// @name Public data members
      
      RunNumber_t    runNumber;     ///< run number
      SubRunNumber_t subRunNumber;  ///< subrun number
      EventNumber_t  eventNumber;   ///< event number
      TimeValue_t    time;          ///< event time
      
      bool           isData;        ///< whether processing real data
      
      std::string processName;      ///< name of the running process
      std::string moduleType;       ///< name of the class of the running module
      std::string moduleLabel;      ///< label of the running module instance
      
      bool isTimeValid;             ///< whether timestamp is valid
      /// @}
      
      
      /// Resets all the fields
      void clear()
        {
          runNumber = 0;
          subRunNumber = 0;
          eventNumber = 0;
          time = 0;
          isData = false;
          processName.clear();
          moduleType.clear();
          moduleLabel.clear();
          isTimeValid = false;
        } // clear()
      
    }; // class EventSeedInputData
  } // namespace SeedServiceHelper
} // namespace sim


#endif // LARSIM_RANDOMUTILS_PROVIDERS_EVENTSEEDINPUTDATA_H