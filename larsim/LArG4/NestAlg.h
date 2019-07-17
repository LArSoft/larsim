#ifndef NESTALG_H
#define NESTALG_H

#include "Geant4/G4Types.hh"
#include "Geant4/G4VParticleChange.hh"
#include <map>

class G4MaterialPropertiesTable;
class G4Step;
class G4Track;

namespace CLHEP { class HepRandomEngine; }

class NestAlg {

 public :

  NestAlg(CLHEP::HepRandomEngine& engine);
  NestAlg(double yieldFactor, CLHEP::HepRandomEngine& engine);

  const G4VParticleChange& CalculateIonizationAndScintillation(G4Track const& aTrack,
							       G4Step  const& aStep);

  void   SetScintillationYieldFactor(double const& yf)     { fYieldFactor = yf;       }
  void   SetScintillationExcitationRatio(double const& er) { fExcitationRatio = er;   }
  int    NumberScintillationPhotons() const                { return fNumScintPhotons; }
  int    NumberIonizationElectrons()  const                { return fNumIonElectrons; }
  double EnergyDeposition()           const                { return fEnergyDep;       }

 private:

  G4double GetGasElectronDriftSpeed(G4double efieldinput,G4double density);
  G4double GetLiquidElectronDriftSpeed(double T, double F, G4bool M, G4int Z);
  G4double CalculateElectronLET ( G4double E, G4int Z );
  G4double UnivScreenFunc ( G4double E, G4double Z, G4double A );
  G4int BinomFluct(G4int N0, G4double prob); //function for doing fluctuations
  void InitMatPropValues ( G4MaterialPropertiesTable* nobleElementMat, int z );

  double             fYieldFactor;     ///< turns scint. on/off
  double 	     fExcitationRatio; ///< N_ex/N_i, the dimensionless ratio of initial
                                       ///< excitons to ions
  int    	     fNumScintPhotons; ///< number of photons produced by the step
  int    	     fNumIonElectrons; ///< number of ionization electrons produced by step
  double 	     fEnergyDep;       ///< energy deposited by the step
  G4VParticleChange  fParticleChange;  ///< pointer to G4VParticleChange
  std::map<int,bool> fElementPropInit; ///< map of noble element z to flag
                                       ///< for whether that element's material
                                       ///< properties table has been initialized
  CLHEP::HepRandomEngine& fEngine;     ///< random engine
};

#endif //NESTALG_H
