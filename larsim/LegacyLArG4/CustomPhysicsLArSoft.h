#ifndef larsim_LegacyLArG4_CustomPhysicsLArSoft_h
#define larsim_LegacyLArG4_CustomPhysicsLArSoft_h

// This file is used to specify template specializations for the
// CustomPhysicsFactory class template, specifically for
// LArSoft-defined physics lists.

#include "larsim/LegacyLArG4/CustomPhysicsFactory.hh"
#include "larsim/LegacyLArG4/OpticalPhysics.hh"
namespace detinfo {
  class DetectorPropertiesData;
}

namespace larg4 {
  template <>
  class CustomPhysicsFactory<OpticalPhysics> : public CustomPhysicsFactoryBase {
  public:
    G4VPhysicsConstructor* Build() const override { return new OpticalPhysics{fDetProp}; }
    CustomPhysicsFactory(detinfo::DetectorPropertiesData const& detProp)
      : CustomPhysicsFactoryBase{"Optical"}, fDetProp{detProp}
    {
      // register self in physics table - note, factory is actually registered
      // in static TheCustomPhysicsTable, not the instance created below
      // which just acts to pass information along
      new CustomPhysicsTable(this);
    }

  private:
    detinfo::DetectorPropertiesData const& fDetProp;
  };
}

#endif
