#ifndef larsim_LegacyLArG4_CustomPhysicsBuiltIns_hh
#define larsim_LegacyLArG4_CustomPhysicsBuiltIns_hh

// The parameters of GEANT4 built in physics modules must be supplied
// to their constructors.  This code initiates CustomPhysicsFactories
// to pass the default QGSP_BERT parameters and add built in physics
// to the CustomPhysicsFactory.  Parameters are passed exactly as they
// were in QGSP_BERT (hence apparent inconsistencies are not the fault
// of the current author!
//
// Any LArSoft-specific lists should be added in CustomPhysicsLArSoft
// (e.g. OpticalPhysics).

#include "Geant4/G4ChargeExchangePhysics.hh"
#include "Geant4/G4DecayPhysics.hh"
#include "Geant4/G4EmExtraPhysics.hh"
#include "Geant4/G4EmLivermorePhysics.hh"
#include "Geant4/G4EmStandardPhysics.hh"
#include "Geant4/G4HadronElasticPhysics.hh"
#include "Geant4/G4HadronElasticPhysicsHP.hh"
#include "Geant4/G4HadronElasticPhysicsPHP.hh"
#include "Geant4/G4HadronPhysicsQGSP_BERT.hh"
#include "Geant4/G4HadronPhysicsQGSP_BERT_HP.hh"
#include "Geant4/G4IonPhysics.hh"
#include "Geant4/G4NeutronTrackingCut.hh"
#include "Geant4/G4StoppingPhysics.hh"

#include "larsim/LegacyLArG4/CustomPhysicsFactory.hh"

namespace larg4 {
  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4EmStandardPhysics>::Build() const
  {
    return new G4EmStandardPhysics(0);
  }

  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4EmExtraPhysics>::Build() const
  {
    return new G4EmExtraPhysics("extra EM");
  }

  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4ChargeExchangePhysics>::Build() const
  {
    return new G4ChargeExchangePhysics(2);
  }

  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4DecayPhysics>::Build() const
  {
    return new G4DecayPhysics("Decay", 0);
  }

  //hadron elastic
  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4HadronElasticPhysics>::Build() const
  {
    return new G4HadronElasticPhysics(0);
  }

  //hadron elastic, neutron HP
  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4HadronElasticPhysicsHP>::Build() const
  {
    return new G4HadronElasticPhysicsHP(0);
  }

  //hadron elastic, particle HP
  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4HadronElasticPhysicsPHP>::Build() const
  {
    return new G4HadronElasticPhysicsPHP(0);
  }

  //hadron inelastic
  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4HadronPhysicsQGSP_BERT>::Build() const
  {
    constexpr G4bool quasiElastic{true};
    return new G4HadronPhysicsQGSP_BERT("hadron", quasiElastic);
  }

  //hadron inelastic, neutron HP
  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4HadronPhysicsQGSP_BERT_HP>::Build() const
  {
    return new G4HadronPhysicsQGSP_BERT_HP("hadronHP", true);
  }

  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4StoppingPhysics>::Build() const
  {
    return new G4StoppingPhysics("stopping");
  }

  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4IonPhysics>::Build() const
  {
    return new G4IonPhysics("ion");
  }

  template <>
  G4VPhysicsConstructor* CustomPhysicsFactory<G4NeutronTrackingCut>::Build() const
  {
    return new G4NeutronTrackingCut("Neutron tracking cut", 0);
  }
}

#endif
