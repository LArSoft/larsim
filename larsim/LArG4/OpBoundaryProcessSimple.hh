/**
 * @file   larsim/LArG4/OpBoundaryProcessSimple.hh
 * @author Ben Jones (bjpjones@mit.edu)
 * @date   March 2010
 * @see    larsim/LArG4/OpBoundaryProcessSimple.cxx
 */


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
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Boundary Process Class Definition
////////////////////////////////////////////////////////////////////////



#ifndef OpBoundaryProcessSimple_h
#define OpBoundaryProcessSimple_h 1

#include "Geant4/G4ForceCondition.hh"
#include "Geant4/G4OpticalPhoton.hh"
#include "Geant4/G4ProcessType.hh"
#include "Geant4/G4String.hh"
#include "Geant4/G4Types.hh"
#include "Geant4/G4VDiscreteProcess.hh"
#include "Geant4/Randomize.hh"

class G4ParticleDefinition;
class G4Step;
class G4Track;
class G4VParticleChange;

/////////////////////
// Class Definition
/////////////////////

namespace larg4 {

  // Possible statuses of a particle after each step.
  enum OpBoundaryProcessSimpleStatus {  Undefined, NotAtBoundary, SimpleAbsorbedNoRefl,
          SimpleAbsorbed, SimpleDiffuse, SimpleSpecular, StepTooSmall, NoRINDEX };

  /**
   * @brief Discrete process for reflection and diffusion at optical interfaces.
   * @see `G4VDiscreteProcess`
   *
   * This class invokes a simplified model of optical reflections at
   * boundaries between different materials.  The relevant reflectivities
   * are ultimately read from `detinfo::LArProperties` via
   * `larg4::MaterialPropertiesLoader`.
   *
   * The required parameters are total reflectance
   * (`detinfo::LArProperties::SurfaceReflectances()`)
   * and ratio of diffuse to specular reflectance
   * (`detinfo::LArProperties::SurfaceReflectanceDiffuseFractions()`). Each
   * photon crossing a boundary with a defined reflectance is randomly either
   * reflected or absorbed and killed according to the supplied probability.
   *
   * Every reflected photon with a defined diffuse reflection fraction
   * is then randomly either diffusely or specularly reflected according
   * to the supplied probability.  All materials with no defined
   * reflectance are assumed to be black and absorb all incident photons.
   *
   * This physics process is loaded in `larg4::OpticalPhysics` physics
   * constructor.
   *
   * This class is based on the `G4OpBoundaryProcess` class in Geant4 and was
   * adapted for LArSoft by Ben Jones, MIT, March 2010.
   */
class OpBoundaryProcessSimple : public G4VDiscreteProcess
  {


  public:


    // Constructors and Destructor

    OpBoundaryProcessSimple(const G4String& processName = "OpBoundary",
          G4ProcessType type = fOptical);

    ~OpBoundaryProcessSimple();



  public:

    G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
    // Returns true -> 'is applicable' only for an optical photon.

    G4double GetMeanFreePath(const G4Track& ,
           G4double ,
           G4ForceCondition* condition);
    // Returns infinity; i. e. the process does not limit the step,
    // but sets the 'Forced' condition for the DoIt to be invoked at
    // every step. However, only at a boundary will any action be
    // taken.

    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
            const G4Step&  aStep);
    // This is the method implementing boundary processes.

    OpBoundaryProcessSimpleStatus GetStatus() const;
    // Returns the current status.


  private:

    G4bool G4BooleanRand(const G4double prob) const;
    // Generate a random bool to decide which process to execute

  private:

    OpBoundaryProcessSimpleStatus fTheStatus;
    G4double fCarTolerance;

    int fVerbosity;

  };


  // Inline methods

  inline
  G4bool OpBoundaryProcessSimple::G4BooleanRand(const G4double prob) const
  {
    /* Returns a random boolean variable with the specified probability */
    return (G4UniformRand() < prob);
  }

  inline
  G4bool OpBoundaryProcessSimple::IsApplicable(const G4ParticleDefinition&
                 aParticleType)
  {
    return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
  }

  inline
  OpBoundaryProcessSimpleStatus OpBoundaryProcessSimple::GetStatus() const
  {
    return fTheStatus;
  }

}

#endif /* OpBoundaryProcessSimple_h */
