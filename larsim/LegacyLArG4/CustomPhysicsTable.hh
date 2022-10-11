////////////////////////////////////////////////////////////////////////
/// \file CustomPhysicsTable.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
//
// The custom physics table keeps track of all compiled physics modules and feeds their
// names and constructors to the ConfigurablePhysicsList.
//
// A CustomPhysicsTable is instantiated in the constructor of each CustomPhysicsFactory<T>
// with a pointer to the instantiating factory as an argument.  This pointer is passed to
// a singleton, static CustomPhysicsTable called TheCustomPhysicsTable.  It is this table
// which the ConfigurablePhysicsList interacts with.
//
// At runtime the CustomPhysicsTable contains an array of CustomPhysicsFactories, one
// for each available physics module, and can provide the list of their names and
// the constructor for each to the physics list.
//
// This means new physics modules can be enabled via the LArG4.xml file with no change to
// the physics list code
//
// See CustomPhysicsFactory.hh and CustomPhysicsFactory.cxx for more details.

#ifndef CUSTOMPHYSICSTABLE_hh
#define CUSTOMPHYSICSTABLE_hh 1

#include <map>
#include <vector>

class G4VPhysicsConstructor;

namespace larg4 {
  class CustomPhysicsFactoryBase;

  class CustomPhysicsTable {
  public:
    CustomPhysicsTable() = default;
    CustomPhysicsTable(CustomPhysicsFactoryBase*);
    std::vector<std::string> GetAvailablePhysicsList();
    bool IsPhysicsAvailable(std::string const&);
    G4VPhysicsConstructor* GetPhysicsConstructor(std::string const& physicsName);
    std::map<std::string, CustomPhysicsFactoryBase*> const& GetFullTable() { return theTable; }
    void AddPhysics(CustomPhysicsFactoryBase*);

  private:
    std::map<std::string, CustomPhysicsFactoryBase*> theTable;
  };

}

#endif
