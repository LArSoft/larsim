////////////////////////////////////////////////////////////////////////
// Class:       ISCalcNESTLAr
// Plugin Type: Algorithm
// File:        ISCalcNESTLAr.cxx
// Description:
// Aug. 30 by Mu Wei
////////////////////////////////////////////////////////////////////////

#include "larsim/IonizationScintillation/ISCalcNESTLAr.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <algorithm>

namespace {
  constexpr double LAr_Z{18};
  constexpr double Density_LAr{1.393};

  constexpr double scint_yield{1.0 / (19.5 * CLHEP::eV)};
  constexpr double resolution_scale{0.107}; // Doke 1976
}

namespace larg4 {

  //----------------------------------------------------------------------------
  ISCalcNESTLAr::ISCalcNESTLAr(CLHEP::HepRandomEngine& Engine)
    : fEngine(Engine), fSCE{lar::providerFrom<spacecharge::SpaceChargeService>()}
  {
    std::cout << "ISCalcNESTLAr Initialize." << std::endl;
  }

  //----------------------------------------------------------------------------
  ISCalcData ISCalcNESTLAr::CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                                            sim::SimEnergyDeposit const& edep)
  {
    CLHEP::RandGauss GaussGen(fEngine);
    CLHEP::RandFlat UniformGen(fEngine);

    double yieldFactor = 1.0; // default quenching factor, for electronic recoils
    double excitationRatio =
      0.21; // ratio for light particle in LAr, such as e-, mu-, Aprile et. al book

    double const energyDeposit = edep.Energy();
    if (energyDeposit < 1 * CLHEP::eV) // too small energy deposition
    {
      return {0., 0., 0., 0.};
    }

    int pdgcode = edep.PdgCode();

    geo::Length_t startx = edep.StartX();
    geo::Length_t starty = edep.StartY();
    geo::Length_t startz = edep.StartZ();
    geo::Length_t endx = edep.EndX();
    geo::Length_t endy = edep.EndY();
    geo::Length_t endz = edep.EndZ();

    double DokeBirks[3];

    double eField = EFieldAtStep(detProp.Efield(), edep);
    if (eField) {
      DokeBirks[0] = 0.07 * pow((eField / 1.0e3), -0.85);
      DokeBirks[2] = 0.00;
    }
    else {
      DokeBirks[0] = 0.0003;
      DokeBirks[2] = 0.75;
    }

    double Density = detProp.Density() /
                     (CLHEP::g / CLHEP::cm3); // argon density at the temperature from Temperature()

    // nuclear recoil quenching "L" factor: total yield is
    // reduced for nuclear recoil as per Lindhard theory
    double epsilon = 11.5 * (energyDeposit / CLHEP::keV) * pow(LAr_Z, (-7. / 3.));

    if (pdgcode == 2112 || pdgcode == -2112) //nuclear recoil
    {
      yieldFactor = 0.23 * (1 + exp(-5 * epsilon)); //liquid argon L_eff
      excitationRatio = 0.69337 + 0.3065 * exp(-0.008806 * pow(eField, 0.76313));
    }

    // determine ultimate number of quanta from current E-deposition (ph+e-) total mean number of exc/ions
    //the total number of either quanta produced is equal to product of the
    //work function, the energy deposited, and yield reduction, for NR
    double MeanNumQuanta = scint_yield * energyDeposit;
    double sigma = sqrt(resolution_scale * MeanNumQuanta); //Fano
    int NumQuanta = int(floor(GaussGen.fire(MeanNumQuanta, sigma) + 0.5));
    double LeffVar = GaussGen.fire(yieldFactor, 0.25 * yieldFactor);
    LeffVar = std::clamp(LeffVar, 0., 1.);

    if (yieldFactor < 1) //nuclear reocils
    {
      NumQuanta = BinomFluct(NumQuanta, LeffVar);
    }

    //if Edep below work function, can't make any quanta, and if NumQuanta
    //less than zero because Gaussian fluctuated low, update to zero
    if (energyDeposit < 1 / scint_yield || NumQuanta < 0) { NumQuanta = 0; }

    // next section binomially assigns quanta to excitons and ions
    int NumExcitons = BinomFluct(NumQuanta, excitationRatio / (1 + excitationRatio));
    int NumIons = NumQuanta - NumExcitons;

    // this section calculates recombination following the modified Birks'Law of Doke, deposition by deposition,
    // may be overridden later in code if a low enough energy necessitates switching to the
    // Thomas-Imel box model for recombination instead (determined by site)
    double dE = energyDeposit / CLHEP::MeV;
    double dx = 0.0;
    double LET = 0.0;
    double recombProb;

    if (pdgcode != 11 && pdgcode != -11 && pdgcode != 13 &&
        pdgcode != -13) //e-: 11, e+: -11, mu-: 13, mu+: -13
    {
      //in other words, if it's a gamma,ion,proton,alpha,pion,et al. do not
      //use the step length provided by Geant4 because it's not relevant,
      //instead calculate an estimated LET and range of the electrons that
      //would have been produced if Geant4 could track them
      LET = CalcElectronLET(1000 * dE);

      if (LET) {
        dx = dE / (Density * LET); //find the range based on the LET
      }

      if (abs(pdgcode) == 2112) //nuclear recoils
      {
        dx = 0;
      }
    }
    else //normal case of an e-/+ energy deposition recorded by Geant
    {
      dx = std::hypot(startx - endx, starty - endy, startz - endz) / CLHEP::cm;
      if (dx) {
        LET = (dE / dx) * (1 / Density); //lin. energy xfer (prop. to dE/dx)
      }
      if (LET > 0 && dE > 0 && dx > 0) {
        double ratio = CalcElectronLET(dE * 1e3) / LET;
        if (ratio < 0.7 && pdgcode == 11) {
          dx /= ratio;
          LET *= ratio;
        }
      }
    }

    DokeBirks[1] = DokeBirks[0] / (1 - DokeBirks[2]); //B=A/(1-C) (see paper)    r
    recombProb = (DokeBirks[0] * LET) / (1 + DokeBirks[1] * LET) +
                 DokeBirks[2]; //Doke/Birks' Law as spelled out in the NEST pape
    recombProb *= (Density / Density_LAr);

    //check against unphysicality resulting from rounding errors
    recombProb = std::clamp(recombProb, 0., 1.);

    //use binomial distribution to assign photons, electrons, where photons
    //are excitons plus recombined ionization electrons, while final
    //collected electrons are the "escape" (non-recombined) electrons
    int const NumPhotons = NumExcitons + BinomFluct(NumIons, recombProb);
    int const NumElectrons = NumQuanta - NumPhotons;

    return {energyDeposit,
            static_cast<double>(NumElectrons),
            static_cast<double>(NumPhotons),
            GetScintYieldRatio(edep)};
  }

  //----------------------------------------------------------------------------
  int ISCalcNESTLAr::BinomFluct(int N0, double prob)
  {
    CLHEP::RandGauss GaussGen(fEngine);
    CLHEP::RandFlat UniformGen(fEngine);

    double mean = N0 * prob;
    double sigma = sqrt(N0 * prob * (1 - prob));
    int N1 = 0;
    if (prob == 0.00) { return N1; }
    if (prob == 1.00) { return N0; }

    if (N0 < 10) {
      for (int i = 0; i < N0; i++) {
        if (UniformGen.fire() < prob) { N1++; }
      }
    }
    else {
      N1 = int(floor(GaussGen.fire(mean, sigma) + 0.5));
    }
    if (N1 > N0) { N1 = N0; }
    if (N1 < 0) { N1 = 0; }
    return N1;
  }

  //----------------------------------------------------------------------------
  double ISCalcNESTLAr::CalcElectronLET(double E)
  {
    double LET;

    if (E >= 1) {
      LET = 116.70 - 162.97 * log10(E) + 99.361 * pow(log10(E), 2) - 33.405 * pow(log10(E), 3) +
            6.5069 * pow(log10(E), 4) - 0.69334 * pow(log10(E), 5) + .031563 * pow(log10(E), 6);
    }
    else if (E > 0 && E < 1) {
      LET = 100;
    }
    else {
      LET = 0;
    }

    return LET;
  }

  //----------------------------------------------------------------------------
  double ISCalcNESTLAr::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
  {
    geo::Point_t pos = edep.MidPoint();
    double EField = efield;
    geo::Vector_t eFieldOffsets;
    if (fSCE->EnableSimEfieldSCE()) {
      eFieldOffsets = fSCE->GetEfieldOffsets(pos);
      EField =
        std::sqrt((efield + efield * eFieldOffsets.X()) * (efield + efield * eFieldOffsets.X()) +
                  (efield * eFieldOffsets.Y() * efield * eFieldOffsets.Y()) +
                  (efield * eFieldOffsets.Z() * efield * eFieldOffsets.Z()));
    }
    return EField;
  }

}
