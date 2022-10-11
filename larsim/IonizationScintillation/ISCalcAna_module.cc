////////////////////////////////////////////////////////////////////////
// Class:       ISCalcAna
// Plugin Type: analyzer (art v2_05_00)
// File:        ISCalcAna_module.cc
//
// Generated at Tue Mar  7 14:59:03 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
// modified at Sept 30, 2019, by muve
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "TNtuple.h"
#include "art_root_io/TFileService.h"

#include "larsim/IonizationScintillation/ISCalc.h"
#include "larsim/IonizationScintillation/ISCalcCorrelated.h"
#include "larsim/IonizationScintillation/ISCalcNESTLAr.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larsim/Simulation/LArG4Parameters.h"

namespace larg4 {
  class ISCalcAna : public art::EDAnalyzer {
  public:
    explicit ISCalcAna(fhicl::ParameterSet const& p);
    ISCalcAna(ISCalcAna const&) = delete;
    ISCalcAna(ISCalcAna&&) = delete;
    ISCalcAna& operator=(ISCalcAna const&) = delete;
    ISCalcAna& operator=(ISCalcAna&&) = delete;

    void analyze(art::Event const& event) override;
    void beginJob() override;
    void endJob() override;

  private:
    std::unique_ptr<ISCalc> fISAlg;
    art::InputTag fEDepTag;
    art::InputTag calcTag; // name of calculator to use, NEST or Separate
    CLHEP::HepRandomEngine& fEngine;
    TNtuple* fNtuple;
  };

  ISCalcAna::ISCalcAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fEDepTag{pset.get<art::InputTag>("SimulationLabel")}
    , calcTag{pset.get<art::InputTag>("ISCalcAlg")}
    , fEngine(art::ServiceHandle<rndm::NuRandomService>()
                ->createEngine(*this, "HepJamesRandom", "NEST", pset, "SeedNEST"))
  {
    std::cout << "ISCalcAna constructor." << std::endl;
  }

  void ISCalcAna::beginJob()
  {
    std::cout << "ISCalcAna beginJob." << std::endl;
    std::cout << "Using " << calcTag.label() << " algorithm to calculate IS." << std::endl;

    if (calcTag.label() == "Separate")
      fISAlg = std::make_unique<larg4::ISCalcSeparate>();
    else if (calcTag.label() == "Correlated") {
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
      fISAlg = std::make_unique<larg4::ISCalcCorrelated>(detProp, fEngine);
    }
    else if (calcTag.label() == "NEST")
      fISAlg = std::make_unique<larg4::ISCalcNESTLAr>(fEngine);
    else
      mf::LogWarning("IonAndScint") << "No ISCalculation set, this can't be good.";

    art::ServiceHandle<art::TFileService const> tfs;
    fNtuple = tfs->make<TNtuple>(
      "nt_is",
      "EDep IS Calc Ntuple",
      "run:event:t:x:y:z:ds:e:trackid:pdg:e_deposit:n_electron:n_photon:scintyield");
  }
  void ISCalcAna::endJob() { std::cout << "ISCalcAna endJob." << std::endl; }

  void ISCalcAna::analyze(art::Event const& event)
  {
    art::Handle<std::vector<sim::SimEnergyDeposit>> edep_handle;
    if (!event.getByLabel(fEDepTag, edep_handle)) {
      std::cout << "PDFastSimPAR Module Cannot getByLabel: " << fEDepTag << std::endl;
      return;
    }

    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);
    for (auto const& edepi : *edep_handle) {
      auto const [energyDeposit, nElectrons, nPhotons, scintYieldRatio] =
        fISAlg->CalcIonAndScint(detProp, edepi);
      fNtuple->Fill(event.run(),
                    event.event(),
                    edepi.T(),
                    edepi.X(),
                    edepi.Y(),
                    edepi.Z(),
                    edepi.StepLength(),
                    edepi.Energy(),
                    edepi.TrackID(),
                    edepi.PdgCode(),
                    energyDeposit,
                    nElectrons,
                    nPhotons,
                    scintYieldRatio);
    }

    std::cout << "ISCalcAna analyze completed." << std::endl;
  }
}
DEFINE_ART_MODULE(larg4::ISCalcAna)
