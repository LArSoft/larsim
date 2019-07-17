////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationNEST.cxx
/// \brief Interface to algorithm class for calculating ionization electrons
///        and scintillation photons using nest
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geant4/G4Step.hh"

#include "larsim/LArG4/ISCalculationNEST.h"
#include "larsim/LArG4/NestAlg.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4{

  //----------------------------------------------------------------------------
  ISCalculationNEST::ISCalculationNEST(CLHEP::HepRandomEngine& engine)
    : fNest(0)
    , fEngine(engine)
  {
    return;
  }

  //----------------------------------------------------------------------------
  ISCalculationNEST::~ISCalculationNEST()
  {
    if(fNest) delete fNest;

    return;
  }

  //----------------------------------------------------------------------------
  void ISCalculationNEST::Initialize()
  {
    // \todo should ideally make the yield factor passed to the NestAlg ctor a parameter
    if(!fNest) fNest = new NestAlg(1., fEngine);

    // Set the step size to small value if NEST is chosen, per Matthew Szydagis,
    // "because without delta rays, the yields are wrong.  The ICARUS model that is
    // in LArSoft uses a fudge factor to compensate, but NEST is "purer" -- no
    // fudge factor. "
    fStepSize = 0.05 * CLHEP::micrometer;

    return;
  }

  //----------------------------------------------------------------------------
  void ISCalculationNEST::Reset()
  {
    fEnergyDeposit   = 0.;
    fNumIonElectrons = 0.;
    fNumScintPhotons = 0.;
    fVisibleEnergyDeposition = 0.;

    return;
  }

  //----------------------------------------------------------------------------
  void ISCalculationNEST::CalculateIonizationAndScintillation(const G4Step* step)
  {
    // get a const representation of the track for this step
    const G4Track track(*(step->GetTrack()));

    fNest->CalculateIonizationAndScintillation(track, *step);

    // compare the energy deposition of this step to what is in the fNest object
    if(fNest->EnergyDeposition() != step->GetTotalEnergyDeposit()/CLHEP::MeV)
      mf::LogWarning("ISCalculationNest") << "NEST and G4 step depositions do not agree!\n"
					  << fNest->EnergyDeposition() << " vs "
					  << step->GetTotalEnergyDeposit()/CLHEP::MeV;

    // Nest uses Geant units, LArSoft assumes energy is in units of MeV here
    fEnergyDeposit   = fNest->EnergyDeposition()/CLHEP::MeV;
    fNumIonElectrons = fNest->NumberIonizationElectrons();
    fNumScintPhotons = fNest->NumberScintillationPhotons();
    fVisibleEnergyDeposition = 0.; //Not implimented for NEST.

    return;
  }

}// namespace
