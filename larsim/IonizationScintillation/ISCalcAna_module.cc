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

#include "art_root_io/TFileService.h"
#include "TNtuple.h"

#include "larsim/IonizationScintillation/ISCalc.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "larsim/IonizationScintillation/ISCalcNESTLAr.h"

#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace larg4 
{
    class ISCalcAna : public art::EDAnalyzer 
    {
    public:
        explicit ISCalcAna(fhicl::ParameterSet const & p);
        ISCalcAna(ISCalcAna const &)                     = delete;
        ISCalcAna(ISCalcAna &&)                          = delete;
        ISCalcAna & operator = (ISCalcAna const &)  = delete;
        ISCalcAna & operator = (ISCalcAna &&)             = delete;
        
        void analyze(art::Event const & event)     override;
        void beginJob()                            override;
        void endJob()                              override;
        
    private:
        ISCalc*                  fISAlg;
        art::InputTag            fEDepTag;
        art::InputTag            calcTag;             // name of calculator to use, NEST or Separate
        CLHEP::HepRandomEngine&  fEngine;
        TNtuple*                 fNtuple;
    };
    
    ISCalcAna::ISCalcAna(fhicl::ParameterSet const & pset)
    : EDAnalyzer(pset)
    , fEDepTag{pset.get<art::InputTag>("SimulationLabel")}
    , calcTag{pset.get<art::InputTag>("ISCalcAlg")}
    , fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "NEST", pset, "SeedNEST"))
    {
        std::cout << "ISCalcAna constructor." << std::endl;
    }
    
    void ISCalcAna::beginJob()
    {
        std::cout << "ISCalcAna beginJob." << std::endl;       
        std::cout << "Using " << calcTag.label() << " algorithm to calculate IS." << std::endl;
        
        if(calcTag.label().compare("NEST") == 0)
        {
            fISAlg = new ISCalcNESTLAr(fEngine);
        }
        else if(calcTag.label().compare("Separate") == 0)
        {
            fISAlg = new ISCalcSeparate();
        }
        else
        {
            mf::LogWarning("ISCalcAna") << "No ISCalculation set, this can't be good.";
        }
        
        fISAlg->Reset();
        fISAlg->Initialize();
        
        art::ServiceHandle<art::TFileService const> tfs;
        fNtuple = tfs->make<TNtuple>("nt_is",
                                     "EDep IS Calc Ntuple",
                                     "run:event:t:x:y:z:ds:e:trackid:pdg:e_deposit:n_electron:n_photon:scintyield");
        return;
    }
    void ISCalcAna::endJob()
    {
        std::cout << "ISCalcAna endJob." << std::endl;
        
        if(fISAlg)
        {
            delete fISAlg;
        }
        
        return;
    }
    void ISCalcAna::analyze(art::Event const & event)
    {
        art::Handle< std::vector<sim::SimEnergyDeposit> > edep_handle;
        if (!event.getByLabel(fEDepTag, edep_handle))
        {
            std::cout << "PDFastSimPAR Module Cannot getByLabel: " << fEDepTag << std::endl;
            return;
        }
        auto const& edeps = edep_handle;
        for(auto const& edepi: *edeps)
        {
            fISAlg->Reset();
            fISAlg->CalcIonAndScint(edepi);
            fNtuple->Fill(event.run(), event.event(), edepi.T(),
                          edepi.X(), edepi.Y(), edepi.Z(), edepi.StepLength(),
                          edepi.Energy(), edepi.TrackID(), edepi.PdgCode(),
                          fISAlg->EnergyDeposit(),
                          fISAlg->NumOfElectrons(),
                          fISAlg->NumOfPhotons(),
//                          fISAlg->NumOfFastPhotons(),
//                          fISAlg->NumOfSlowPhotons());
                          fISAlg->ScintillationYieldRatio());
        }
        
        std::cout << "ISCalcAna analyze completed." << std::endl;
        return;
    }
}
DEFINE_ART_MODULE(larg4::ISCalcAna)
