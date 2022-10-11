////////////////////////////////////////////////////////////////////////
/// \file CustomPhysicsFactory.h
//
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
#ifndef CUSTOMPHYSICSFACTORY_hh
#define CUSTOMPHYSICSFACTORY_hh

#include "Geant4/G4VPhysicsConstructor.hh"
#include "larsim/LegacyLArG4/CustomPhysicsTable.hh"

namespace larg4 {
  class CustomPhysicsTable;

  class CustomPhysicsFactoryBase {
  public:
    explicit CustomPhysicsFactoryBase(std::string const& name) : fName{name}
    {
      if (empty(name))
        std::cerr << "CustomPhysicsFactory Error : Physics registered with no name!" << std::endl;
    }
    virtual ~CustomPhysicsFactoryBase() = default;

    std::string const& GetName() const noexcept { return fName; }

    virtual G4VPhysicsConstructor* Build() const = 0;

  private:
    std::string fName;
  };

  template <class T>
  class CustomPhysicsFactory : public CustomPhysicsFactoryBase {
  public:
    explicit CustomPhysicsFactory(std::string const& name);

    G4VPhysicsConstructor* Build() const { return new T{}; }
  };

  template <class T>
  CustomPhysicsFactory<T>::CustomPhysicsFactory(std::string const& name)
    : CustomPhysicsFactoryBase{name}
  {
    // register self in physics table - note, factory is actually registered
    // in static TheCustomPhysicsTable, not the instance created below
    // which just acts to pass information along
    new CustomPhysicsTable(this);
  }
}

#endif
