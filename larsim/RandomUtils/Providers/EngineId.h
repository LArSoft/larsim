/**
 * @file EngineId.h
 * @brief An identifier for random engines
 * @author Rob Kutschke (kutschke@fnal.gov)
 * 
 * An identifier may consist of simply a module label or a module label plus an
 * instance name.
 */

#ifndef LARSIM_RANDOMUTILS_PROVIDERS_ENGINEID_H
#define LARSIM_RANDOMUTILS_PROVIDERS_ENGINEID_H 1

#include <string>
#include <ostream>

namespace sim {

  /// Namespace for implementation details of SeedMaster
  namespace SeedMasterHelper {
    
    /// Identifier for a engine, made of module name and optional instance name
    struct EngineId {
      
      /// structure to identify a "global" flavour constructor
      struct Global_t {};
      
      /// A constant to select a "global" flavour constructor
      static Global_t global;
      
      /// Constructor (module name is required)
      EngineId(std::string const& mod, std::string const& inst = std::string()):
        moduleLabel(mod),
        instanceName(inst)
        {}
      
      /// Constructor (module name is required)
      EngineId(std::string const& inst, Global_t):
        moduleLabel(),
        instanceName(inst)
        {}
      
      // Accept compiler written d'tor, copy c'tor, copy and move assignments.
      
      /// Returns whether the label is "global" (no module context)
      bool isGlobal() const { return moduleLabel.empty(); }
      
      /// Returns whether the instance label is defined
      bool hasInstanceName() const { return !instanceName.empty(); }
      
      /// Sets this ID to the specified global instance
      void setGlobal(std::string inst)
        { moduleLabel.clear(); instanceName = inst; }
      
      
      /// Returns true if both module and instance names match
      bool operator== (EngineId const& rhs) const
        {
          if ( moduleLabel  != rhs.moduleLabel  ) return false;
          if ( instanceName != rhs.instanceName ) return false;
          return true;
        } // operator== ()
      
      /// Lexicographic sort (module name first, then instance name)
      bool operator< (EngineId const& rhs) const
        {
          if (moduleLabel  < rhs.moduleLabel) return true;
          if (moduleLabel == rhs.moduleLabel) {
            if (instanceName < rhs.instanceName) return true;
          }
          return false;
        } // operator< ()
      
      /// Converts the information in a module_name[.instance_name] string
      operator std::string() const
        {
          std::string id = moduleLabel;
          if (hasInstanceName()) id.append(1, '.').append(instanceName);
          return id;
        } // operator std::string()
      
      /// Converts the information in a module_name:instance_name string
      std::string artName() const { return moduleLabel + ':' + instanceName; }
      
      std::string moduleLabel; ///< module label
      std::string instanceName; ///< instance name
      
    }; // end class EngineId
    
    inline std::ostream& operator<<
      (std::ostream& ost, const EngineId& id )
      { return ost << std::string(id); }
    
  } // end namespace SeedMasterHelper
  
} // namespace sim

#endif // LARSIM_RANDOMUTILS_PROVIDERS_ENGINEID_H
