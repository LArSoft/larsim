////////////////////////////////////////////////////////////////////////
/// \file   AuxDetReadout.h
/// \brief  A Geant4 sensitive detector that accumulates information.
/// \author miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_AUXDETREADOUT_H
#define LARG4_AUXDETREADOUT_H

#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/globals.hh"

#include "Simulation/AuxDetSimChannel.h"
#include "Geometry/Geometry.h"
#include "Geometry/AuxDetGeo.h"
#include "Geometry/AuxDetSensitiveGeo.h"

#include <vector>

// Forward declarations
class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

namespace larg4 {
  
  class AuxDetReadout : public G4VSensitiveDetector
  {
  public:
    // Constructor.
    AuxDetReadout(std::string const& name, 
		  unsigned int       adNum,
		  unsigned int       svNum);
    
    // Destructor
    virtual ~AuxDetReadout();
    
    // Required for classes that inherit from G4VSensitiveDetector.
    //
    // Called at start and end of each event.
    virtual void Initialize(G4HCofThisEvent*);
    virtual void EndOfEvent(G4HCofThisEvent*);
    
    // Called to clear any accumulated information.
    virtual void clear();
    
    // The key method of this class.  It's called by Geant4 for each
    // step within the read-out geometry.  It accumulates the energy
    // in the G4Step in the ?.
    virtual G4bool ProcessHits( G4Step*, G4TouchableHistory* );
    
    // Moved here from AuxDetSimChannel
    virtual void AddParticleStep(int	inputTrackID,
				 float	inputEnergyDeposited,
				 float	inputEntryX,
				 float	inputEntryY,
				 float	inputEntryZ,
				 float	inputEntryT,
				 float	inputExitX,
				 float	inputExitY,
				 float	inputExitZ,
				 float	inputExitT,
				 float	inputExitMomentumX,
				 float	inputExitMomentumY,
				 float	inputExitMomentumZ);

    // Empty methods; they have to be defined, but they're rarely
    // used in Geant4 applications.
    virtual void DrawAll();
    virtual void PrintAll();
    
    // Independent method; returns the accumulated information
    sim::AuxDetSimChannel const GetAuxDetSimChannel() const { return fAuxDetSimChannel; };
    
  private:
    art::ServiceHandle<geo::Geometry> fGeoHandle;        ///< Handle to the Geometry service
    uint32_t                          fAuxDet;           ///< which AuxDet this AuxDetReadout corresponds to
    uint32_t                          fAuxDetSensitive;  ///< which sensitive volume of the AuxDet this AuxDetReadout corresponds to
    sim::AuxDetSimChannel             fAuxDetSimChannel; ///< Contains the sim::AuxDetSimChannel for this AuxDet
    std::vector<sim::AuxDetIDE>       fAuxDetIDEs;       ///< list of IDEs in one channel
};
}

#endif // LARG4_AUXDETREADOUT_H
