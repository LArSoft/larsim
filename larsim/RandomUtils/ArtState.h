/**
 * @file ArtState.h
 * @brief Class holding the current state of art processing
 * @author Rob Kutschke (kutschke@fnal.gov)
 */

#ifndef LARSIM_RANDOMUTILS_ARTSTATE_H
#define LARSIM_RANDOMUTILS_ARTSTATE_H 1

// LArSeedService libraries
#include "larsim/RandomUtils/Providers/EventSeedInputData.h"

// framework libraries
#include "art/Utilities/Exception.h"
#include "art/Persistency/Provenance/EventAuxiliary.h"
#include "art/Persistency/Provenance/Timestamp.h"
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Framework/Principal/Event.h"

// supporting libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <string>


namespace sim {
  
  namespace SeedServiceHelper {
    
    /// Describe the current state of art processing, as understood by the LArSeedService.
    class ArtState {
        public:
      typedef enum {
        unDefined,             ///< not assigned yet
        inServiceConstructor,  ///< in service construction phase
        inModuleConstructor,   ///< in module construction phase
        inBeginRun,            ///< in begin of run phase
        inModuleBeginRun,      ///< in begin of run for a module
        inEvent,               ///< in event phase
        inModuleEvent,         ///< in event processing by a module
        inEndJob,              ///< in end job
        inOther                ///< none of the above
      } state_type; ///< type of state of art (what's doing)
      
      using EventInfo_t = art::EventAuxiliary;
      
      
      ArtState(state_type start_state = unDefined)
        : artState(start_state)
        , lastEvent()
        , lastModule()
        , procName()
        {}
      
      // Accept compiler written d'tor, copy c'tor and copy assignment.
      
      //@{
      /// Records the status of ART
      void set_state(state_type astate) { artState = astate; }
      
      /// Records the new status of ART and returns the old one
      state_type transit_to(state_type astate)
        {
          state_type old_state = state();
          set_state(astate);
          LOG_DEBUG("ArtState")
            << "LArSeedService::ArtState: transition from "
            << stateName(old_state) << " to " << stateName();
          return old_state;
        } // transit_to()
      
      /// Resets the status to "something else" (inOther)
      void reset_state() { transit_to(inOther); }
      
      /// Records the specified event ID
      void set_event(art::Event const& evt)
        {
          lastEvent = {
            evt.id(), evt.time(), evt.isRealData(), evt.experimentType()
            };
        } // set_event()
      void reset_event() { lastEvent = EventInfo_t(); }
      
      /// Records the specified module description
      void set_module(art::ModuleDescription const& desc)
        {
          lastModule = desc;
          set_process_name(desc);
        }
      void reset_module() { lastModule = art::ModuleDescription(); }
      
      void set_process_name(std::string pn) { procName = pn; }
      void set_process_name(art::ModuleDescription const& currentModuleDesc)
        {
          if (!procName.empty()
            && (procName != currentModuleDesc.processName()))
          {
            throw art::Exception(art::errors::LogicError)
              << "Process name changed from '" << procName << "' to '"
              << currentModuleDesc.processName() << "'?!?";
          }
          procName = currentModuleDesc.processName();
        } // set_process_name()
      
      //@}
      
      
      //@{
      /// Getters
      
      state_type state() const { return artState; }
      
      std::string stateName() const { return stateName(state()); }
      
      art::EventID const& eventID() const { return lastEvent.id(); }
      
      EventInfo_t const& eventInfo() const { return lastEvent; }
      
      art::ModuleDescription const& moduleDesc() const { return lastModule; }
      
      std::string moduleLabel() const { return lastModule.moduleLabel(); }
      
      std::string processName() const { return procName; }
      //@}
      
      
      EventSeedInputData getEventSeedInputData() const
        {
          EventSeedInputData data;
          data.runNumber = eventID().run();
          data.subRunNumber = eventID().subRun();
          data.eventNumber = eventID().event();
          
          data.time = eventInfo().time().value();
          data.isTimeValid
            = (eventInfo().time() != art::Timestamp::invalidTimestamp());
          
          data.isData = eventInfo().isRealData();
          
          data.processName = processName();
          data.moduleType = moduleDesc().moduleName();
          data.moduleLabel = moduleLabel();
          
          return data;
        } // getEventSeedInputData()
      /// @}
      
      static std::string stateName(state_type state)
        {
          switch (state) {
            case unDefined:            return "(not assigned yet)";
            case inServiceConstructor: return "service construction";
            case inModuleConstructor:  return "module construction";
            case inBeginRun:           return "begin of run";
            case inModuleBeginRun:     return "begin of run for module";
            case inEvent:              return "event preparation";
            case inModuleEvent:        return "event processing by a module";
            case inEndJob:             return "end job";
            case inOther:              return "unidentified";
          } // switch
          throw art::Exception(art::errors::LogicError)
            << "artext::SeedServiceHelper::ArtState::stateName: unknown state #"
            << ((int) state) << "\n";
        } // stateName()
      
        protected:
      state_type artState; ///< current state of the art
      
      EventInfo_t lastEvent;
      art::ModuleDescription lastModule;
      std::string procName;
    }; // end ArtState

  } // end namespace SeedServiceHelper

} // namespace sim

#endif // LARSIM_RANDOMUTILS_ARTSTATE_H
