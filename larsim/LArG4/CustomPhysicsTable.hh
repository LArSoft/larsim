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

#ifndef CUSTOMPHYSICSFACTORY_hh
namespace larg4 { class CustomPhysicsFactoryBase; }
#endif

#include "Geant4/G4VPhysicsConstructor.hh"


namespace larg4 {
  class CustomPhysicsTable
  {
    
  public:
    CustomPhysicsTable(CustomPhysicsFactoryBase*);
    ~CustomPhysicsTable() {};
    std::vector<std::string> GetAvailablePhysicsList();
    bool IsPhysicsAvailable(std::string);
    G4VPhysicsConstructor * GetPhysicsConstructor(std::string);
    std::map<std::string,CustomPhysicsFactoryBase* > GetFullTable()
    { return theTable;}
    void AddPhysics(CustomPhysicsFactoryBase*);
    
    
  protected:
    std::map<std::string,CustomPhysicsFactoryBase* > theTable;
    CustomPhysicsTable() {};
  };
  
}

#include "larsim/LArG4/CustomPhysicsFactory.hh"
  
#endif




// Sept 2009 - Ben Jones, MIT
