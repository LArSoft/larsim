////////////////////////////////////////////////////////////////////////
/// \file  ISCalculation.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larsim/LegacyLArG4/ISCalculation.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "Geant4/G4Step.hh"

namespace larg4 {

  //......................................................................
  double ISCalculation::EFieldAtStep(double efield, const G4Step* step) const
  {
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    if (!SCE->EnableSimEfieldSCE()) return efield;
    geo::Point_t midPoint{
      (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition()) * 0.5 /
      CLHEP::cm};
    auto const eFieldOffsets = SCE->GetEfieldOffsets(midPoint);
    return efield * std::hypot(1 + eFieldOffsets.X(), eFieldOffsets.Y(), eFieldOffsets.Z());
  }

}
