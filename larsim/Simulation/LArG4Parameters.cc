////////////////////////////////////////////////////////////////////////
/// \file  LArG4Parameters_service.cc
/// \brief Store parameters for running LArG4
///
/// \author bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// This service exists to pass parameters to various different
// classes in LArG4, which are not necessary directly called by
// the LArG4_module class.
//
// Ben Jones, MIT, March 2010

#include "larsim/Simulation/LArG4Parameters.h"

#include "fhiclcpp/ParameterSet.h"

namespace sim {

  //--------------------------------------------------------------------------
  LArG4Parameters::LArG4Parameters(fhicl::ParameterSet const& pset)
    : fOpVerbosity{pset.get<int>("OpticalSimVerbosity")}
    , fParticleKineticECut{pset.get<double>("ParticleKineticEnergyCut")}
    , fStoreTrajectories{pset.get<bool>("StoreTrajectories")}
    , fDrawNeutrals{pset.get<bool>("VisualizeNeutrals")}
    , fVisualizationEnergyCut{pset.get<double>("VisualizationEnergyCut")}
    , fUseCustomPhysics{pset.get<bool>("UseCustomPhysics")}
    , fModifyProtonCut{pset.get<bool>("ModifyProtonCut", false)}
    , fNewProtonCut{fModifyProtonCut ? pset.get<double>("NewProtonCut") /* for HadronHP */ : 0.0}
    , fLongitudinalDiffusion{pset.get<double>("LongitudinalDiffusion")}
    , fTransverseDiffusion{pset.get<double>("TransverseDiffusion")}
    , fElectronClusterSize{pset.get<double>("ElectronClusterSize")}
    , fMinNumberOfElCluster{pset.get<int>("MinNumberOfElCluster")}
    , fEnabledPhysics{pset.get<std::vector<std::string>>("EnabledPhysics")}
    , fK0Bias{pset.get<int>("CosmogenicK0Bias")}
    , fXSBias{pset.get<int>("CosmogenicXSMNBiasFactor")}
    , fXBias{pset.get<int>("CosmogenicXSMNBiasOn")}
    // First of last 3 flags above turns on secondary particle bias for
    // K0s,Lambdas,neutrons in MuNuclear.
    // The 2nd is the enhancement factor for XS bias in MuNuclear. Keep it
    // <=100.
    // The 3rd turns on cross-section bias in MuNuclear.
    , fKeepEMShowerDaughters{pset.get<bool>("KeepEMShowerDaughters")}
    , fDisableWireplanes{pset.get<bool>("DisableWireplanes")}
    , fSkipWireSignalInTPCs{pset.get<std::vector<unsigned short int>>("SkipWireSignalInTPCs")}
    , fRecombA{pset.get<double>("RecombA", util::kRecombA)}
    , fRecombk{pset.get<double>("Recombk", util::kRecombk)}
    , fModBoxA{pset.get<double>("ModBoxA", util::kModBoxA)}
    , fModBoxB{pset.get<double>("ModBoxB", util::kModBoxB)}
    , fLarqlChi0A{pset.get<double>("LarqlChi0A")}
    , fLarqlChi0B{pset.get<double>("LarqlChi0B")}
    , fLarqlChi0C{pset.get<double>("LarqlChi0C")}
    , fLarqlChi0D{pset.get<double>("LarqlChi0D")}
    , fLarqlAlpha{pset.get<double>("LarqlAlpha")}
    , fLarqlBeta{pset.get<double>("LarqlBeta")}
    , fWph{pset.get<double>("Wph")}
    , fUseModBoxRecomb{pset.get<bool>("UseModBoxRecomb")}
    , fUseModLarqlRecomb{pset.get<bool>("UseModLarqlRecomb")}
    , fUseBinomialFlucts{pset.get<bool>("UseBinomialFlucts", true)}
    , fIonAndScintCalculator{pset.get<std::string>("IonAndScintCalculator", "Separate")}
    , fOpticalParamVolumes{pset.get<std::vector<std::string>>("OpticalParamVolumes")}
    , fOpticalParamModels{pset.get<std::vector<std::string>>("OpticalParamModels")}
    , fOpticalParamOrientations{pset.get<std::vector<int>>("OpticalParamOrientations")}
    , fOpticalParamParameters{pset.get<std::vector<std::vector<std::vector<double>>>>(
        "OpticalParamParameters")}
    , fLitePhotons{pset.get<bool>("UseLitePhotons")}
    , fFillSimEnergyDeposits{pset.get<bool>("FillSimEnergyDeposits", false)}
    , fNoElectronPropagation{pset.get<bool>("NoElectronPropagation", false)}
    , fNoPhotonPropagation{pset.get<bool>("NoPhotonPropagation", false)}
  {}
}
