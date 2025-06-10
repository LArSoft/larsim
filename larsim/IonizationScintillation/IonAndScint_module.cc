////////////////////////////////////////////////////////////////////////
// Class:       IonAndScint
// Plugin Type: producer
// File:        IonAndScint_module.cc
// Description:
// - acts on sim::SimEnergyDeposit from LArG4Main,
// - calculate the number of photons and electrons
// Input: 'sim::SimEnergyDeposit'
// Output: updated 'sim::SimEnergyDeposit' with numPhotons and numElectrons
//
//This module calculate the number of photons and electrons produced at each step where energy is deposited.
//The Separate algorithm is used by default, but this can be changed via the "ISCalcAlg"
//fhicl parameter tag.
//At the end of this module the numPhotons and numElectrons of sim:SimEnergyDeposit have been updated.
//
// Aug.18 by Mu Wei
//
// 10/28/2019 Wenqiang Gu (wgu@bnl.gov)
//            Add the Space Charge Effect (SCE) if the option is enabled
////////////////////////////////////////////////////////////////////////

// LArSoft includes

#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/IonizationScintillation/ISCalc.h"
#include "larsim/IonizationScintillation/ISCalcCorrelated.h"
#include "larsim/IonizationScintillation/ISCalcNESTLAr.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "larsim/IonizationScintillation/IGapChargeTransport.hh"

#include "nurandom/RandomUtils/NuRandomService.h"

// Framework includes

#include "art/Utilities/make_tool.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandomEngine.h"

#include <iostream>
#include <sstream> // std::stringstream
#include <string>
#include <vector>

using std::string;
using SimEnergyDepositCollection = std::vector<sim::SimEnergyDeposit>;

//make_tool function for a tool moving the charge from the CRP gap (DUNE only)
namespace {
  std::unique_ptr<gap::IGapChargeTransport>
  make_gap_tool(fhicl::ParameterSet const& tool_pset) {
  if (tool_pset.is_empty()) return nullptr;
   std::cout << "Building the charge recovery tool!" << std::endl; 
   return art::make_tool<gap::IGapChargeTransport>(tool_pset);
  }
}

namespace larg4 {
  class IonAndScint : public art::EDProducer {
  public:
    explicit IonAndScint(fhicl::ParameterSet const& pset);
    void produce(art::Event& event) override;
    void beginJob() override;
    void endJob() override;

  private:
    std::vector<art::Handle<SimEnergyDepositCollection>> inputCollections(art::Event const&) const;

    // name of calculator: Separate, Correlated, or NEST
    art::InputTag calcTag;

    // The input module labels to specify which SimEnergyDeposit
    // collections to use as input. Specify only the module label,
    // not the instance name(s). Instances to be used have to be
    // passed via the Instances parameter.
    // If empty, this module uses all the available collections.
    std::vector<std::string> fInputModuleLabels;

    std::unique_ptr<ISCalc> fISAlg;
    CLHEP::HepRandomEngine& fEngine;
    string Instances;
    std::vector<string> instanceNames;
    bool fSavePriorSCE;
    std::unique_ptr<gap::IGapChargeTransport> fGapTool;
  };

  //......................................................................
  IonAndScint::IonAndScint(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , calcTag{pset.get<art::InputTag>("ISCalcAlg")}
    , fInputModuleLabels{pset.get<std::vector<std::string>>("InputModuleLabels", {})}
    , fEngine(art::ServiceHandle<rndm::NuRandomService> {}
              -> registerAndSeedEngine(createEngine(0, "HepJamesRandom", "ISCalcAlg"),
                                       "HepJamesRandom",
                                       "ISCalcAlg",
                                       pset,
                                       "SeedISCalcAlg"))
    , Instances{pset.get<string>("Instances", "LArG4DetectorServicevolTPCActive")}
    , fSavePriorSCE{pset.get<bool>("SavePriorSCE", false)}
    , fGapTool{make_gap_tool(pset.get<fhicl::ParameterSet>("ChargeTransport", {}))}
  {
    std::cout << "IonAndScint Module Construct" << std::endl;
    std::cout << "==== ParameterSet Keys for IonAndScint ====" << std::endl;
    for (const auto& name : pset.get_names()) {
      std::cout << " - " << name << std::endl;
    }
    if (Instances.empty()) {
      std::cout << "Produce SimEnergyDeposit in default volume - LArG4DetectorServicevolTPCActive"
                << std::endl;
      instanceNames.push_back("LArG4DetectorServicevolTPCActive");
    }
    else {
      std::stringstream input(Instances);
      string temp;
      while (std::getline(input, temp, ';')) {
        instanceNames.push_back(temp);
      }

      std::cout << "Produce SimEnergyDeposit in volumes: " << std::endl;
      for (auto instanceName : instanceNames) {
        std::cout << " - " << instanceName << std::endl;
      }
    }

    produces<std::vector<sim::SimEnergyDeposit>>();
    if (fSavePriorSCE) produces<std::vector<sim::SimEnergyDeposit>>("priorSCE");

    if (fGapTool) {
      std::cout << "fGapTool is configured" << std::endl;
    } else {
      //std::cout << "fGapTool is nullptr (tool not configured)" << std::endl;
    }
  }

  //......................................................................
  void IonAndScint::beginJob()
  {
    std::cout << "IonAndScint beginJob." << std::endl;
    std::cout << "Using " << calcTag.label() << " algorithm to calculate IS." << std::endl;

    if (calcTag.label() == "Separate")
      fISAlg = std::make_unique<ISCalcSeparate>();
    else if (calcTag.label() == "Correlated") {
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
      fISAlg = std::make_unique<ISCalcCorrelated>(detProp, fEngine);
    }
    else if (calcTag.label() == "NEST")
      fISAlg = std::make_unique<ISCalcNESTLAr>(fEngine);
    else
      mf::LogWarning("IonAndScint") << "No ISCalculation set, this can't be good.";
    if (fGapTool) std::cout << "Tool on, begin job" << std::endl; 
  }

  //......................................................................
  void IonAndScint::endJob()
  {
    std::cout << "IonAndScint endJob." << std::endl;
  }

  //......................................................................
  std::vector<art::Handle<SimEnergyDepositCollection>> IonAndScint::inputCollections(
    art::Event const& e) const
  {
    if (empty(fInputModuleLabels)) {
      mf::LogDebug("IonAndScint") << "Retrieving all products" << std::endl;
      return e.getMany<SimEnergyDepositCollection>();
    }

    std::vector<art::Handle<SimEnergyDepositCollection>> result;

    for (auto const& module : fInputModuleLabels) {

      mf::LogDebug("IonAndScint") << "Retrieving products with module label " << module
                                  << std::endl;

      auto handels = e.getMany<SimEnergyDepositCollection>(art::ModuleLabelSelector(module));

      if (empty(handels)) {
        throw art::Exception(art::errors::ProductNotFound)
          << "IonAndScint module cannot find any SimEnergyDeposits with module label " << module
          << " as requested in InputModuleLabels. \n";
      }

      result.insert(result.end(), handels.begin(), handels.end());
    }

    return result;
  }

  //......................................................................
  void IonAndScint::produce(art::Event& event)
  {
    std::cout << "IonAndScint Module Producer" << std::endl;

    std::vector<art::Handle<SimEnergyDepositCollection>> edepHandle = inputCollections(event);

    if (empty(edepHandle)) {
      std::cout << "IonAndScint Module Cannot Retrive SimEnergyDeposit" << std::endl;
      return;
    }

    auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);

    auto simedep = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
    auto simedep1 = std::make_unique<std::vector<sim::SimEnergyDeposit>>(); // for prior-SCE depos
    for (auto edeps : edepHandle) {
      // Do some checking before we proceed
      if (!edeps.isValid()) {
        std::cout << "!edeps.isValid()" << std::endl;
        continue;
      }

      auto index = std::find(
        instanceNames.begin(), instanceNames.end(), edeps.provenance()->productInstanceName());
      if (index == instanceNames.end()) {
        std::cout << "Skip SimEnergyDeposit in: " << edeps.provenance()->productInstanceName()
                  << std::endl;
        continue;
      }

      std::cout << "SimEnergyDeposit input module: " << edeps.provenance()->moduleLabel()
                << ", instance name: " << edeps.provenance()->productInstanceName() << std::endl;

      for (sim::SimEnergyDeposit const& edepi : *edeps) {
        auto const isCalcData = fISAlg->CalcIonAndScint(detProp, edepi);

        int ph_num = round(isCalcData.numPhotons);
        int ion_num = round(isCalcData.numElectrons);
        float scintyield = isCalcData.scintillationYieldRatio;
        float edep_tmp = edepi.Energy();
        geo::Point_t startPos_tmp = edepi.Start();
        geo::Point_t endPos_tmp = edepi.End();
        double startTime_tmp = edepi.StartT();
        double endTime_tmp = edepi.EndT();
        int trackID_tmp = edepi.TrackID();
        int pdgCode_tmp = edepi.PdgCode();
        int origTrackID_tmp = edepi.OrigTrackID();
	int ion_num_tmp = ion_num;
        if (sce->EnableSimSpatialSCE()) {
          auto posOffsetsStart =
            sce->GetPosOffsets({edepi.StartX(), edepi.StartY(), edepi.StartZ()});
          auto posOffsetsEnd = sce->GetPosOffsets({edepi.EndX(), edepi.EndY(), edepi.EndZ()});
          startPos_tmp =
            geo::Point_t{(float)(edepi.StartX() - posOffsetsStart.X()), //x should be subtracted
                         (float)(edepi.StartY() + posOffsetsStart.Y()),
                         (float)(edepi.StartZ() + posOffsetsStart.Z())};
          endPos_tmp =
            geo::Point_t{(float)(edepi.EndX() - posOffsetsEnd.X()), //x should be subtracted
                         (float)(edepi.EndY() + posOffsetsEnd.Y()),
                         (float)(edepi.EndZ() + posOffsetsEnd.Z())};
        }


       
	//----------gap charge-----------------------
	if (fGapTool){
	if (fGapTool->EnableGapChargeTransport()){
           if (edeps.provenance()->productInstanceName() == "LArG4DetectorServicevolEnclosureTPC"){
               std::pair<geo::Point_t, int> pair_ = fGapTool->GetOffset(edepi.StartX(), edepi.StartY(), edepi.StartZ(), ion_num);
	       ion_num_tmp = pair_.second;
	       startPos_tmp = pair_.first; 
               endPos_tmp = geo::Point_t{edepi.EndX() + startPos_tmp.X() - edepi.StartX(), 
			         edepi.EndY() + startPos_tmp.Y() - edepi.StartY(),
				 edepi.EndZ() + startPos_tmp.Z() - edepi.StartZ()};

		std::cout << "Cherge moved to: " << startPos_tmp.X() << "; " << startPos_tmp.Y()  << "; " << startPos_tmp.Z() << std::endl;
           }
        }
        }
	
        if (!fGapTool || (fGapTool &&!(fGapTool->EnableGapChargeTransport() && ion_num_tmp <= 0))) simedep->emplace_back(ph_num,
                              ion_num_tmp,
                              scintyield,
                              edep_tmp,
                              startPos_tmp,
                              endPos_tmp,
                              startTime_tmp,
                              endTime_tmp,
                              trackID_tmp,
                              pdgCode_tmp,
                              origTrackID_tmp);
	//----------gap charge-----------------------
	       

	if (fSavePriorSCE) {
          simedep1->emplace_back(ph_num,
                                 ion_num,
                                 scintyield,
                                 edepi.Energy(),
                                 edepi.Start(),
                                 edepi.End(),
                                 edepi.StartT(),
                                 edepi.EndT(),
                                 edepi.TrackID(),
                                 edepi.PdgCode(),
                                 edepi.OrigTrackID());
        }
      }
    }
    event.put(std::move(simedep));
    if (fSavePriorSCE) event.put(std::move(simedep1), "priorSCE");
  }
} // namespace
DEFINE_ART_MODULE(larg4::IonAndScint)
