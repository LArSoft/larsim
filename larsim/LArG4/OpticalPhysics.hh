// Modified for LArG4 by Ben Jones, MIT, Sept 2009
//
// This class is based upon the default Geant4 optical physics constructor.
// However, for optical physics in LArSoft:
//
// - Cerenkov light can be enabled or disabled as required
//
// - The default geant4 reflection process is swapped for a simpler
//     custom boundary process in LArG4
//
// Parameters relating to the physics in this constructor can be set
// in the Utilities/LArParameters service.
//
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
// GEANT4 tag $Name:  $
//
//---------------------------------------------------------------------------
//
// ClassName:   OpticalPhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 05.12.2005 V.Ivanchenko add controlled verbosity
// 23.11.2006 V.Ivanchenko remove mscStepLimit option and improve cout
//
//----------------------------------------------------------------------------
//
// This class provides construction of default EM standard physics
//

#ifndef OpticalPhysics_h
#define OpticalPhysics_h 1

#include "Geant4/G4String.hh"
#include "Geant4/G4Types.hh"
#include "Geant4/G4VPhysicsConstructor.hh"
class G4Cerenkov;
class G4OpAbsorption;
class G4OpBoundaryProcess;
class G4OpRayleigh;
class G4OpWLS;
class G4Scintillation;

#if defined __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-private-field"
#endif

namespace larg4 {

  class OpBoundaryProcessSimple;

  class OpticalPhysics : public G4VPhysicsConstructor
  {
  public:
    OpticalPhysics(G4int ver = 0, const G4String& name = "G4Optical");
    virtual ~OpticalPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4int  verbose;
    G4Cerenkov*              fTheCerenkovProcess;
    G4Scintillation*         fTheScintillationProcess;
    G4OpAbsorption*          fTheAbsorptionProcess;
    G4OpRayleigh*            fTheRayleighScatteringProcess;
    OpBoundaryProcessSimple* fTheBoundaryProcess;
    G4OpBoundaryProcess*     fTheBoundaryProcess_g4;
    G4OpWLS *                fTheWLSProcess;
  };

}
#if defined __clang__
  #pragma clang diagnostic pop
#endif

#endif
