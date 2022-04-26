////////////////////////////////////////////////////////////////////////
// Class:       ISCalcCorrelated
// Plugin Type: algorithm
// File:        ISCalcCorrelated.h and ISCalcCorrelated.cxx
// Description: Interface to algorithm class for a specific calculation of
//              ionization electrons and scintillation photons, based on
//              simple microphysics arguments to establish an anticorrelation
//              between these two quantities.
// Input: 'sim::SimEnergyDeposit'
// Output: num of Photons and Electrons
// May 2020 by W Foreman
// Modified: Adding corrections for low electric field (LArQL model)
// Jun 2020 by L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////

#ifndef IS_ISCALCCORRELATED_H
#define IS_ISCALCCORRELATED_H

#include "larsim/IonizationScintillation/ISCalc.h"
#include "larsim/IonizationScintillation/ISTPC.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandEngine.h"


namespace larg4 {
  class ISCalcCorrelated : public ISCalc {
  public:
    ISCalcCorrelated(detinfo::DetectorPropertiesData const& detProp, CLHEP::HepRandomEngine& Engine);

    double EFieldAtStep(double efield,
                        sim::SimEnergyDeposit const& edep)
      override; //value of field with any corrections for this step
    ISCalcData CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                               sim::SimEnergyDeposit const& edep) override;

  private:
    double fGeVToElectrons;   ///< from LArG4Parameters service
    double fWion;             ///< W_ion (23.6 eV) == 1/fGeVToElectrons
    double fWph;              ///< from LArG4Parameters service
    double fScintPreScale;    ///< scintillation pre-scaling factor from LArProperties service
    double fRecombA;          ///< from LArG4Parameters service
    double fRecombk;          ///< from LArG4Parameters service
    double fModBoxA;          ///< from LArG4Parameters service
    double fModBoxB;          ///< from LArG4Parameters service
    double fLarqlChi0A;       ///< from LArG4Parameters service
    double fLarqlChi0B;       ///< from LArG4Parameters service
    double fLarqlChi0C;       ///< from LArG4Parameters service
    double fLarqlChi0D;       ///< from LArG4Parameters service
    double fLarqlAlpha;       ///< from LArG4Parameters service
    double fLarqlBeta;        ///< from LArG4Parameters service
    bool fUseModBoxRecomb;    ///< from LArG4Parameters service
    bool fUseModLarqlRecomb;  ///< from LArG4Parameters service

    const spacecharge::SpaceCharge* fSCE;
    const detinfo::LArProperties* fLArProp;

    void CalcIon(sim::SimEnergyDeposit const& edep);
    void CalcScint(sim::SimEnergyDeposit const& edep);
    double EscapingEFraction(double const dEdx); //LArQL chi0 function = fraction of escaping electrons
    double FieldCorrection(double const EF, double const dEdx); //LArQL f_corr function = correction factor for electric field dependence
    int BinomFluct(int N0, double prob);

    ISTPC fISTPC;
    CLHEP::HepRandomEngine& fEngine; // random engine 
  };
}
#endif
