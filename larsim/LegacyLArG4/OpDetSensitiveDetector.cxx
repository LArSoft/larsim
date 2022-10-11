////////////////////////////////////////////////////////////////////////
/// \file OpDetSensitiveDetector.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Implementation of the OpDetSensitiveDetector
//
// See comments in OpDetSensitiveDetector.h
//
// Ben Jones, MIT, 06/04/2010
//

#include "larsim/LegacyLArG4/OpDetSensitiveDetector.h"
#include "Geant4/G4SDManager.hh"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/LegacyLArG4/OpDetLookup.h"
#include "larsim/LegacyLArG4/OpDetPhotonTable.h"

namespace {

  /// Converts a photon `energy` [eV] into its Wavelength [nm]
  constexpr double Wavelength(double energy);

} // local namespace

namespace larg4 {

  OpDetSensitiveDetector::OpDetSensitiveDetector(G4String DetectorUniqueName,
                                                 bool useLitePhotons /* = false */)
    : G4VSensitiveDetector(DetectorUniqueName), fUseLitePhotons(useLitePhotons)
  {
    // Register self with sensitive detector manager
    G4SDManager::GetSDMpointer()->AddNewDetector(this);

    // Get instances of singleton classes
    fTheOpDetLookup = OpDetLookup::Instance();
    fThePhotonTable = OpDetPhotonTable::Instance();
  }

  //--------------------------------------------------------

  void OpDetSensitiveDetector::AddLitePhoton(G4Step const* aStep, int OpDet)
  {

    double const time = aStep->GetTrack()->GetGlobalTime();

    // the guideline: if it's VUV (~128 nm) is direct, otherwise it is reflected
    double const energy = aStep->GetTrack()->GetVertexKineticEnergy() / CLHEP::eV;
    bool const reflected = Wavelength(energy) > 200.0; // nm

    // Add this photon to the detected photons table
    fThePhotonTable->AddLitePhoton(OpDet, static_cast<int>(time), 1, reflected);

  } // OpDetSensitiveDetector::AddLitePhoton()

  //--------------------------------------------------------

  void OpDetSensitiveDetector::AddPhoton(G4Step const* aStep, int OpDet)
  {
    sim::OnePhoton ThePhoton;

    // Get photon data to store in the hit

    ThePhoton.SetInSD = true;

    auto const& track = *(aStep->GetTrack());
    auto const& startPos = track.GetVertexPosition();
    ThePhoton.InitialPosition = {startPos.x(), startPos.y(), startPos.z()};

    //ThePhoton.Time                = track.GetGlobalTime() - fGlobalTimeOffset;
    ThePhoton.Time = track.GetGlobalTime();

    ThePhoton.Energy = track.GetVertexKineticEnergy();

    // Lookup which OpDet we are in
    G4StepPoint const* preStepPoint = aStep->GetPreStepPoint();

    // Store relative position on the photon detector
    G4ThreeVector worldPosition = preStepPoint->GetPosition();
    G4ThreeVector localPosition =
      preStepPoint->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(
        worldPosition);
    ThePhoton.FinalLocalPosition = {
      localPosition.x() / CLHEP::cm, localPosition.y() / CLHEP::cm, localPosition.z() / CLHEP::cm};

    // Add this photon to the detected photons table
    fThePhotonTable->AddPhoton(OpDet, std::move(ThePhoton));

  } // OpDetSensitiveDetector::AddPhoton()

  //--------------------------------------------------------

  G4bool OpDetSensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
  {
    // Lookup which OpDet we are in
    int const OpDet = fTheOpDetLookup->GetOpDet(aStep->GetPreStepPoint()->GetPhysicalVolume());

    // Add this photon to the detected photons table
    if (fUseLitePhotons)
      AddLitePhoton(aStep, OpDet);
    else
      AddPhoton(aStep, OpDet);

    // Kill this photon track
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);

    return true;
  }

  //--------------------------------------------------------

  void OpDetSensitiveDetector::Initialize(G4HCofThisEvent*) {}

}

//--------------------------------------------------------
namespace {

  constexpr double Wavelength(double energy)
  {

    // SI 2019 (eV nm):
    constexpr double hc = 6.62607015e-34 * 299792458.0 / 1.602176634e-19 * 1e9;

    return hc / energy; // nm

  } // Wavelength()

} // local namespace

//--------------------------------------------------------
