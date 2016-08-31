////////////////////////////////////////////////////////////////////////
/// \file  ISCalculation.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larsim/LArG4/ISCalculation.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// Framework includes                                                                                                         
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace larg4 {

  //----------------------------------------------------------------------                                                    
  ISCalculation::ISCalculation()
  {
  }

  //----------------------------------------------------------------------                                                    
  ISCalculation::~ISCalculation()
  {
  }

  //......................................................................                                                    
  double ISCalculation::EFieldAtStep(double fEfield, const G4Step* step) const
  {
    double EField = fEfield;
    std::vector<double> EfieldOffsets;
    CLHEP::Hep3Vector EfieldVec;
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    if (SCE->EnableSimEfieldSCE())
      {
        G4ThreeVector midPoint = 0.5*( step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition() );
        EfieldOffsets = SCE->GetEfieldOffsets(midPoint.x()/CLHEP::cm,midPoint.y()/CLHEP::cm,midPoint.z()/CLHEP::cm);
        EfieldVec.set(fEfield + fEfield*EfieldOffsets.at(0), fEfield*EfieldOffsets.at(1), fEfield*EfieldOffsets.at(2));
        EField = EfieldVec.mag();
      }

    return EField;
  }

}
