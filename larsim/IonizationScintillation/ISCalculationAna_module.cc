////////////////////////////////////////////////////////////////////////
// Class:       ISCalculationAna
// Plugin Type: analyzer (art v2_05_00)
// File:        ISCalculationAna_module.cc
//
// Generated at Tue Mar  7 14:59:03 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
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

#include "art/Framework/Services/Optional/TFileService.h"
#include "TNtuple.h"

#include "ISCalculationSeparate.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

namespace larg4 {
  class ISCalculationAna;
}


class larg4::ISCalculationAna : public art::EDAnalyzer {
public:
  explicit ISCalculationAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ISCalculationAna(ISCalculationAna const &) = delete;
  ISCalculationAna(ISCalculationAna &&) = delete;
  ISCalculationAna & operator = (ISCalculationAna const &) = delete;
  ISCalculationAna & operator = (ISCalculationAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p);

private:

  ISCalculationSeparate fISAlg;
  art::InputTag         fEDepTag;

  TNtuple* fNtuple;
};


larg4::ISCalculationAna::ISCalculationAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  this->reconfigure(p);
}

void larg4::ISCalculationAna::analyze(art::Event const & e)
{
  auto const& edep_handle = e.getValidHandle< std::vector<sim::SimEnergyDeposit> >(fEDepTag);
  auto const& edep_vec(*edep_handle);

  for(auto const& edep : edep_vec){
    fISAlg.Reset();
    fISAlg.CalculateIonizationAndScintillation(edep);
    fNtuple->Fill(e.run(),e.event(),
		  edep.T(),
		  edep.X(),edep.Y(),edep.Z(),edep.StepLength(),
		  edep.Energy(),edep.TrackID(),edep.PdgCode(),
		  fISAlg.EnergyDeposit(),
		  fISAlg.NumberIonizationElectrons(),
		  fISAlg.NumberScintillationPhotons());
  }
}

void larg4::ISCalculationAna::beginJob()
{
  art::ServiceHandle<sim::LArG4Parameters> lgpHandle;
  fISAlg.Initialize(lar::providerFrom<detinfo::LArPropertiesService>(),
		    lar::providerFrom<detinfo::DetectorPropertiesService>(),
		    &(*lgpHandle),
		    lar::providerFrom<spacecharge::SpaceChargeService>());
  art::ServiceHandle<art::TFileService> tfs;
  fNtuple = tfs->make<TNtuple>("nt_is","EDep IS Calc Ntuple","run:event:t:x:y:z:ds:e:trackid:pdg:e_deposit:n_electron:n_photon");
}

void larg4::ISCalculationAna::reconfigure(fhicl::ParameterSet const & p)
{
  fEDepTag = p.get<art::InputTag>("EDepModuleLabel");
}

DEFINE_ART_MODULE(larg4::ISCalculationAna)
