////////////////////////////////////////////////////////////////////////
/// \file OpParamSD.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
//
// This class represents a partially opaque, parameterized optical volume.
//
// The photon transmission probability is determined as a function of
// photon position and momentum direction by a derived class of 
// OpParamAction. The particlar derived class to use is specified in a string given to 
// the constructor, along with any parameters required to define that object.
// This implementation allows for generic extensions to new types of wireplanes
// or opaque surfaces for future liquid argon TPC detectors. 
//
// On each step of a photon within this volume, the GetAttenuationFraction
// method of the OpParamAction derivative is called to get a photon 
// transmission probability.  A fraction of all photons are killed 
// accordingly.
//
// This sensitive detector object is attached to physical volumes by the 
// OpDetReadoutGeometry class.
//
//
// Ben Jones, MIT,2013
//


// LArSoft includes

#include "Geant4/G4VSensitiveDetector.hh"
#include "larsim/Simulation/sim.h"
#include "larsim/LArG4/OpParamAction.h"


// Geant4 includes

#include "Geant4/Randomize.hh"
#include "Geant4/G4RandomTools.hh"


#ifndef OpParamSD_h
#define OpParamSD_h 1

class G4HCofThisEvent;
class G4TOuchableHistory;
class G4Step;

namespace sim{
  class SimPhotonsCollection;
}

namespace larg4 {

  class OpDetLookup;
  class OpDetPhotonTable;

  class OpParamSD : public G4VSensitiveDetector
  {

    
  public:
    OpParamSD(G4String name, std::string ModelName, int Orientation, std::vector<std::vector<double> > Parameters);
    virtual ~OpParamSD(){}
    
    
    // Beginning and end of event
    virtual void Initialize(G4HCofThisEvent*);
    virtual void EndOfEvent(G4HCofThisEvent*){}
    
    // Tidy up event in abort
    virtual void clear(){}
    
    // Run per step in sensitive volume
    virtual G4bool ProcessHits( G4Step*, G4TouchableHistory*);
    
    
    // Required but empty
    virtual void DrawAll(){}
    virtual void PrintAll(){}

    
  private:
    G4bool G4BooleanRand(const G4double prob) const;
    
    OpParamAction *        fOpa;
    std::map<G4int, bool>  fPhotonAlreadyCrossed; 
    
  };



  inline
    G4bool OpParamSD::G4BooleanRand(const G4double prob) const
    {
      /* Returns a random boolean variable with the specified probability */
      return (G4UniformRand() < prob);
    }



}

#endif
