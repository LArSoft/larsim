////////////////////////////////////////////////////////////////////////
/// \file  LArG4Parameters.h
/// \brief Store parameters for running LArG4
///
/// \author bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
//
// This service exists to pass parameters to various different
// classes in LArG4, which are not necessary directly called by
// the LArG4_module class.
//
// Ben Jones, MIT, March 2010

#ifndef LArG4Parameters_h
#define LArG4Parameters_h 1

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

namespace fhicl {
  class ParameterSet;
}

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

#include <string>
#include <vector>

namespace sim {

  class LArG4Parameters {
  public:
    LArG4Parameters(fhicl::ParameterSet const& pset);

    int OpVerbosity() const { return fOpVerbosity; }
    double ParticleKineticEnergyCut() const { return fParticleKineticECut; }
    bool StoreTrajectories() const { return fStoreTrajectories; }
    bool DrawNeutrals() const { return fDrawNeutrals; }
    double VisualizationEnergyCut() const { return fVisualizationEnergyCut; }
    bool UseCustomPhysics() const { return fUseCustomPhysics; }
    bool ModifyProtonCut() const { return fModifyProtonCut; }
    double NewProtonCut() const { return fNewProtonCut; }
    double RecombA() const { return fRecombA; }
    double Recombk() const { return fRecombk; }
    double ModBoxA() const { return fModBoxA; }
    double ModBoxB() const { return fModBoxB; }
    double EllipsModBoxA() const { return fEllipsModBoxA; }
    double EllipsModBoxB() const { return fEllipsModBoxB; }
    double EllipsModBoxR() const { return fEllipsModBoxR; }
    double LarqlChi0A() const { return fLarqlChi0A; }
    double LarqlChi0B() const { return fLarqlChi0B; }
    double LarqlChi0C() const { return fLarqlChi0C; }
    double LarqlChi0D() const { return fLarqlChi0D; }
    double LarqlAlpha() const { return fLarqlAlpha; }
    double LarqlBeta() const { return fLarqlBeta; }
    double Wph() const { return fWph; }
    double QAlpha() const { return fQAlpha; }
    bool UseModBoxRecomb() const { return fUseModBoxRecomb; }
    bool UseEllipsModBoxRecomb() const { return fUseEllipsModBoxRecomb; }
    bool UseModLarqlRecomb() const { return fUseModLarqlRecomb; }
    bool UseBinomialFlucts() const { return fUseBinomialFlucts; }
    double GeVToElectrons() const { return util::kGeVToElectrons; }
    double LongitudinalDiffusion() const { return fLongitudinalDiffusion; }
    double TransverseDiffusion() const { return fTransverseDiffusion; }
    double ElectronClusterSize() const { return fElectronClusterSize; }
    int MinNumberOfElCluster() const { return fMinNumberOfElCluster; }
    const std::vector<std::string>& EnabledPhysics() const { return fEnabledPhysics; }
    int K0Bias() const { return fK0Bias; }
    int MNXBias() const { return fXBias; }
    int MNXSBias() const { return fXSBias; }
    bool KeepEMShowerDaughters() const { return fKeepEMShowerDaughters; }
    bool DisableWireplanes() const { return fDisableWireplanes; }
    const std::vector<unsigned short int>& SkipWireSignalInTPCs() const
    {
      return fSkipWireSignalInTPCs;
    }
    const std::string& IonAndScintCalculator() const { return fIonAndScintCalculator; }
    const std::vector<std::string>& OpticalParamVolumes() const { return fOpticalParamVolumes; }
    const std::vector<std::string>& OpticalParamModels() const { return fOpticalParamModels; }
    const std::vector<int>& OpticalParamOrientations() const { return fOpticalParamOrientations; }
    const std::vector<std::vector<std::vector<double>>>& OpticalParamParameters() const
    {
      return fOpticalParamParameters;
    }
    bool UseLitePhotons() const { return fLitePhotons; }

    bool FillSimEnergyDeposits() const { return fFillSimEnergyDeposits; }
    bool NoElectronPropagation() const { return fNoElectronPropagation; }
    bool NoPhotonPropagation() const { return fNoPhotonPropagation; }

  private:
    int const fOpVerbosity; ///< Verbosity of optical simulation - soon to be depricated
    double const
      fParticleKineticECut; ///< Minimum energy a particle needs in order to be stored in the particle list [GeV]
    bool const fStoreTrajectories;        ///< Whether to store full trajectories for every particle
                                          ///< simulated by Geant4
    bool const fDrawNeutrals;             ///< depricated
    double const fVisualizationEnergyCut; ///< depricated, GeV
    bool const fUseCustomPhysics;         ///< Whether to use a custom list of physics processes
                                          ///< or the default
    bool const fModifyProtonCut; ///< Whether to enable custom ProtonCut value, needed for HadronHP
    double const fNewProtonCut;  ///< New Proton Cut parameter to override default in HadronHP
    double const
      fLongitudinalDiffusion; ///< Amount of diffusion in the longitudinal direction, cm^2/ns
    double const fTransverseDiffusion; ///< Amount of diffusion in the transverse direction, cm^2/ns
    double const fElectronClusterSize; ///< Number of ionization electrons in a given cluster
                                       ///< to be simulated in the readout simulation
    int const fMinNumberOfElCluster;   ///< Minimum number of electron clusters
    std::vector<std::string> const
      fEnabledPhysics; ///< List of enabled physics processes if using Custom physics
    int const fK0Bias; ///< Turns on secondary particle bias for K0, Lambda,
                       ///< neutrons in MuNuclear
    int const fXSBias; ///< Turns on cross-section bian in MuNuclear
    int const fXBias;  ///< Enhancement factor for cross-section bian in MuNuclear,
                       ///< should be <= 100
    bool const fKeepEMShowerDaughters; ///< Whether to keep the secondary, tertiary, etc.
                                       ///< particles from an EM shower in the output
    bool const fDisableWireplanes;     ///< Turn of LAr sensitivity and remove charge
                                       ///< drift simulation - use for running pure optical sims
    std::vector<unsigned short int> const
      fSkipWireSignalInTPCs;     ///< selective disabling of drift simulation
    double const fRecombA;       ///< Possibly override the RecombA parameter
    double const fRecombk;       ///< Possibly override the Recombk parameter
    double const fModBoxA;       ///< Possibly override the ModBoxA parameter
    double const fModBoxB;       ///< Possibly override the ModBoxB parameter
    double const fEllipsModBoxA; ///< Possibly override the EllipsModBoxA parameter
    double const fEllipsModBoxB; ///< Possibly override the EllipsModBoxB parameter
    double const fEllipsModBoxR; ///< Possibly override the EllipsModBoxR parameter
    double const fLarqlChi0A;    ///< Possibly override the LarqlChi0A parameter
    double const fLarqlChi0B;    ///< Possibly override the LarqlChi0B parameter
    double const fLarqlChi0C;    ///< Possibly override the LarqlChi0C parameter
    double const fLarqlChi0D;    ///< Possibly override the LarqlChi0D parameter
    double const fLarqlAlpha;    ///< Possibly override the LarqlAlpha parameter
    double const fLarqlBeta;     ///< Possibly override the LarqlBeta parameter
    double const fWph;           ///< Possibly override the Wph parameter
    double const fQAlpha;        ///< Possibly override the QAlpha parameter
    bool const fUseModBoxRecomb; ///< Use Modified Box model recombination instead of Birks
    bool const
      fUseEllipsModBoxRecomb; ///< Use Ellipsoid Modified Box model recombination instead of Birks
    bool const fUseModLarqlRecomb; ///< Use LArQL model recombination correction (dependence on EF)
    bool const fUseBinomialFlucts; ///< Use binomial fluctuations in correlated method
    std::string const
      fIonAndScintCalculator; ///< Name of algorithm to use to calculate the number of
                              ///< ionization electrons and scintillation photons
                              ///< for each G4 step, used by
                              ///< LArG4/IonizationAndScintillation.cxx
    std::vector<std::string> const
      fOpticalParamVolumes; ///< List of volume names which have parameterized
                            ///< optical models
    std::vector<std::string> const fOpticalParamModels; ///< List of names of those models
    std::vector<int> const
      fOpticalParamOrientations; ///< List of orientations of (eg wireplane) in each
                                 ///< param volume
    std::vector<std::vector<std::vector<double>>> const
      fOpticalParamParameters; ///< Model dependent list of
                               ///< parameters for optically
                               ///< parameterized volumes

    bool const fLitePhotons;

    bool const fFillSimEnergyDeposits; ///< handle to fill SimEdeps or not
    bool const fNoElectronPropagation; ///< specifically prevents electron propagation
    bool const fNoPhotonPropagation;   ///< specifically prevents photon propagation in opfast
  };
}

DECLARE_ART_SERVICE(sim::LArG4Parameters, SHARED)
#endif
