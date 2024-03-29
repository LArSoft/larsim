//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// GEANT4 tag $Name: geant4-09-02 $
//
//
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "larsim/LegacyLArG4/G4ThermalElectron.hh"
#include "Geant4/G4ParticleTable.hh"

// ######################################################################
// ###                    THERMAL (DRIFT) ELECTRON                    ###
// ######################################################################
G4ThermalElectron* G4ThermalElectron::theInstance = 0;

G4ThermalElectron* G4ThermalElectron::Definition()
{
  if (theInstance != 0) return theInstance;
  const G4String name = "thermalelectron";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0) {
    // create particle
    //
    //    Arguments for constructor are as follows
    //               name             mass          width         charge
    //             2*spin           parity  C-conjugation
    //          2*Isospin       2*Isospin3       G-parity
    //               type    lepton number  baryon number   PDG encoding
    //             stable         lifetime    decay table
    //             shortlived      subType    anti_encoding

    // use constants in CLHEP
    //  static const double CLHEP::electron_mass_c2 = 0.51099906 * MeV;

    anInstance = new G4ParticleDefinition(name,
                                          CLHEP::electron_mass_c2,
                                          0.0 * CLHEP::MeV,
                                          -1. * CLHEP::eplus,
                                          1,
                                          0,
                                          0,
                                          0,
                                          0,
                                          0,
                                          "lepton",
                                          1,
                                          0,
                                          11,
                                          true,
                                          -1.0,
                                          NULL,
                                          false,
                                          "e");
    // Bohr Magnetron
    G4double muB =
      -0.5 * CLHEP::eplus * CLHEP::hbar_Planck / (CLHEP::electron_mass_c2 / CLHEP::c_squared);

    anInstance->SetPDGMagneticMoment(muB * 2. * 1.0011596521859);
  }
  theInstance = reinterpret_cast<G4ThermalElectron*>(anInstance);
  return theInstance;
}

G4ThermalElectron* G4ThermalElectron::ThermalElectronDefinition()
{
  return Definition();
}

G4ThermalElectron* G4ThermalElectron::ThermalElectron()
{
  return Definition();
}
