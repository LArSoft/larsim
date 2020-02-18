// Modified for LArG4 by Ben Jones, MIT, Nov 2012
//
// Thisclass is based upon thedefaultGeant4 optical physics constructor.
// However, forfast optical physics inLArSoft:
//
// - The default geant4 scintillation process is swapped for the LArG4 fast
//     scintillation process
//
// - Cerenkov light (still slow) can beenabledor disabled as required
//
// - The default geant4reflection process is swapped for a simpler
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
// ClassName:   FastOpticalPhysics
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

#ifndef FastOpticalPhysics_h
#define FastOpticalPhysics_h 1


#include "Geant4/G4VPhysicsConstructor.hh"
#include "Geant4/G4String.hh"
#include "Geant4/G4Types.hh"

#if defined __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-private-field"
#endif

class G4Cerenkov;
class G4OpAbsorption;
class G4OpWLS;
class G4OpRayleigh;

namespace larg4 {

  class OpFastScintillation;
  class OpBoundaryProcessSimple;

  class FastOpticalPhysics : public G4VPhysicsConstructor
  {
  public:
    FastOpticalPhysics(G4int ver = 0, const G4String& name = "FastOptical");
    virtual ~FastOpticalPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4int  verbose;
    G4Cerenkov*              fTheCerenkovProcess;
    OpFastScintillation*     fTheScintillationProcess;
    G4OpAbsorption*          fTheAbsorptionProcess;
    G4OpRayleigh*            fTheRayleighScatteringProcess;
    OpBoundaryProcessSimple* fTheBoundaryProcess;
    G4OpWLS *                fTheWLSProcess;
  };

}
#if defined __clang__
  #pragma clang diagnostic pop
#endif

#endif
