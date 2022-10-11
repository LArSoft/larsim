//////////////////////////////////////////////////////////////////////////////
/// \file NeutronHPphysics.hh
/// \brief High precision neutron physics constructor for Geant4
/// \details Adapted from the Geant4 example Hadr04
//////////////////////////////////////////////////////////////////////////////
#pragma once

// Geant4 includes
#include "Geant4/G4VPhysicsConstructor.hh"

class NeutronHPphysics : public G4VPhysicsConstructor {
public:
  NeutronHPphysics(const G4String& name = "neutron");

  void ConstructParticle() override {}
  void ConstructProcess() override;

  void SetThermalPhysics(G4bool flag) { fThermal = flag; };

private:
  G4bool fThermal;
};
