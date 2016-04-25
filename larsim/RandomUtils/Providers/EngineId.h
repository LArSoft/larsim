/**
 * @file EngineId.h
 * @brief An identifier for random engines
 * @author Rob Kutschke (kutschke@fnal.gov)
 * 
 * An identifier may consist of simply a module label or a module label plus an
 * instance name.
 */

#ifndef SeedService_EngineId_h
#define SeedService_EngineId_h

#include <string>
#include <ostream>

namespace sim {

  /// Namespace for implementation details of SeedMaster
  namespace SeedMasterHelper {
    
    /// Identifier for a engine, made of module name and optional instance name
    struct EngineId {
      
      /// Constructor (module name is required)
      EngineId(std::string const& mod, std::string const& inst = std::string()):
        moduleLabel(mod),
        instanceName(inst)
        {}
      
      // Accept compiler written d'tor, copy c'tor, copy and move assignments.
      
      /// Returns whether the instance label is defined
      bool hasInstanceName() const { return !instanceName.empty(); }
      
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

#endif /* SeedService_EngineId_h */
