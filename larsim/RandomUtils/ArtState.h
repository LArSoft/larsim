/**
 * @file ArtState.h
 * @brief Class holding the current state of art processing
 * @author Rob Kutschke (kutschke@fnal.gov)
 */

#ifndef SeedService_ArtState_h
#define SeedService_ArtState_h

// C/C++ standard libraries
#include <string>

// framework libraries
#include "art/Utilities/Exception.h"
#include "art/Persistency/Provenance/EventAuxiliary.h"
#include "art/Persistency/Provenance/Timestamp.h"
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Framework/Principal/Event.h"

// LArSeedService libraries
#include "larsim/RandomUtils/Providers/EventSeedInputData.h"


namespace sim {
  
  namespace SeedServiceHelper {
    
    /// Describe the current state of art processing, as understood by the LArSeedService.
    class ArtState {
        public:
      typedef enum {
        unDefined,             ///< not assigned yet
        inModuleConstructor,   ///< in module construction phase
        inBeginRun,            ///< in begin of run phase
        inEvent,               ///< in event phase
        inModuleEvent,         ///< in event processing by a module
        inOther                ///< none of the above
      } state_type; ///< type of state of art (what's doing)
      
      using EventInfo_t = art::EventAuxiliary;
      
      
      ArtState()
        : artState(unDefined)
        , lastEvent()
        , lastModule()
        , procName()
        {}
      
      // Accept compiler written d'tor, copy c'tor and copy assignment.
      
      //@{
      /// Records the status of ART
      void set_state(state_type astate) { artState = astate; }
      
      /// Resets the status to "something else" (inOther)
      void reset_state() { set_state(inOther); }
      
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
      
        protected:
      state_type artState; ///< current state of the art
      
      EventInfo_t lastEvent;
      art::ModuleDescription lastModule;
      std::string procName;
    }; // end ArtState

  } // end namespace SeedServiceHelper

} // namespace sim

#endif /* SeedService_ArtState_h */
