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
#include "larsim/IonizationScintillation/ISCalc.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "larsim/IonizationScintillation/ISCalcNESTLAr.h"
#include "larsim/IonizationScintillation/ISCalcCorrelated.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"

#include <stdio.h>
#include <sstream>      // std::stringstream, std::stringbuf
#include <iostream>
#include <string>

using std::string;

namespace larg4
{
    class IonAndScint : public art::EDProducer
    {
    public:
        explicit IonAndScint(fhicl::ParameterSet const& pset);
        void produce(art::Event& event)         override;
        void beginJob()                         override;
        void endJob()                           override;

    private:
        art::InputTag            calcTag;             // name of calculator to use, NEST or Separate
        ISCalc*                  fISAlg;
        CLHEP::HepRandomEngine&  fEngine;             // random number engine for NEST algorithm
        string                   Instances;
        std::vector<string>      instanceNames;
    };

    //......................................................................
    IonAndScint::IonAndScint(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , calcTag{pset.get<art::InputTag>("ISCalcAlg")}
    , fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "NEST", pset, "SeedNEST"))
    , Instances{pset.get<string>("Instances", "LArG4DetectorServicevolTPCActive"),}
    {
        std::cout << "IonAndScint Module Construct" << std::endl;

        if(Instances.empty())
        {
            std::cout << "Produce SimEnergyDeposit in default volume - LArG4DetectorServicevolTPCActive" << std::endl;
            instanceNames.push_back("LArG4DetectorServicevolTPCActive");
            // produces< std::vector<sim::SimEnergyDeposit> >("LArG4DetectorServicevolTPCActive");
        }
        else
        {
            std::stringstream input(Instances);
            string temp;
            while(std::getline(input, temp, ';'))
            {
                instanceNames.push_back(temp);
            }

            std::cout << "Produce SimEnergyDeposit in volumes: " << std::endl;
            for (auto instanceName : instanceNames)
            {
               std::cout << " - " << instanceName << std::endl;
               // produces< std::vector<sim::SimEnergyDeposit> >(instanceName);
            }
        }

        produces< std::vector<sim::SimEnergyDeposit> >();
    }

    //......................................................................
    void IonAndScint::beginJob()
    {
        std::cout << "IonAndScint beginJob." << std::endl;
        std::cout << "Using " << calcTag.label() << " algorithm to calculate IS." << std::endl;
        if(calcTag.label().compare("NEST") == 0)
        {
            fISAlg = new ISCalcNESTLAr(fEngine);
        }
        else if(calcTag.label().compare("Separate") == 0)
        {
            fISAlg = new ISCalcSeparate();
        }
        else if(calcTag.label().compare("Correlated") == 0)
        {
            fISAlg = new ISCalcCorrelated();
        }
        else
        {
            mf::LogWarning("ISCalcAna") << "No ISCalculation set, this can't be good.";
        }

        fISAlg->Initialize();
        fISAlg->Reset();

        return;
    }

    //......................................................................
    void IonAndScint::endJob()
    {
        std::cout << "IonAndScint endJob." << std::endl;

        if(fISAlg) delete fISAlg;

        return;
    }

    //......................................................................
    void IonAndScint::produce(art::Event& event)
    {
        std::cout << "IonAndScint Module Producer" << std::endl;
        std::vector< art::Handle< std::vector< sim::SimEnergyDeposit > > > edepHandle;

        event.getManyByType(edepHandle);

        if (edepHandle.size() == 0)
        {
            std::cout << "IonAndScint Module Cannot Retrive SimEnergyDeposit" << std::endl;
            return;
        }

        auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();

        std::unique_ptr< std::vector<sim::SimEnergyDeposit> >  simedep (new std::vector<sim::SimEnergyDeposit>);
        for (auto edeps: edepHandle)
        {
            // Do some checking before we proceed
            if (!edeps.isValid())
            {
                std::cout << "!edeps.isValid()" << std::endl;
                continue;
            }

            auto index = std::find(instanceNames.begin(), instanceNames.end(), edeps.provenance()->productInstanceName());
            if (index == instanceNames.end())
            {
                std::cout << "Skip SimEnergyDeposit in: " << edeps.provenance()->productInstanceName() << std::endl;
                continue;
            }

            std::cout << "SimEnergyDeposit input module: " << edeps.provenance()->moduleLabel()
                      << ", instance name: " << edeps.provenance()->productInstanceName() << std::endl;

            // std::unique_ptr< std::vector<sim::SimEnergyDeposit> >  simedep (new std::vector<sim::SimEnergyDeposit>);
            for (sim::SimEnergyDeposit const& edepi: *edeps)
            {
                fISAlg->CalcIonAndScint(edepi);

                int           ph_num        = round(fISAlg->NumOfPhotons());
                int           ion_num       = round(fISAlg->NumOfElectrons());
                float         scintyield    = fISAlg->ScintillationYieldRatio();
                float         edep_tmp      = edepi.Energy();
                geo::Point_t  startPos_tmp  = edepi.Start();
                geo::Point_t  endPos_tmp    = edepi.End();
                double        startTime_tmp = edepi.StartT();
                double        endTime_tmp   = edepi.EndT();
                int           trackID_tmp   = edepi.TrackID();
                int           pdgCode_tmp   = edepi.PdgCode();

                if(sce->EnableSimSpatialSCE()) {
                    auto posOffsetsStart = sce->GetPosOffsets({edepi.StartX(),edepi.StartY(),edepi.StartZ()});
                    auto posOffsetsEnd = sce->GetPosOffsets({edepi.EndX(),edepi.EndY(),edepi.EndZ()});
                    startPos_tmp = geo::Point_t{(float)(edepi.StartX()-posOffsetsStart.X()), //x should be subtracted
                                   (float)(edepi.StartY()+posOffsetsStart.Y()),
                                   (float)(edepi.StartZ()+posOffsetsStart.Z())};
                    endPos_tmp   = geo::Point_t{(float)(edepi.EndX()-posOffsetsEnd.X()), //x should be subtracted
                                   (float)(edepi.EndY()+posOffsetsEnd.Y()),
                                   (float)(edepi.EndZ()+posOffsetsEnd.Z())};
                }

                simedep->emplace_back(ph_num, ion_num, scintyield, edep_tmp,
                                      startPos_tmp, endPos_tmp,
                                      startTime_tmp, endTime_tmp,
                                      trackID_tmp, pdgCode_tmp);
            }

            // event.put(std::move(simedep), edeps.provenance()->productInstanceName().c_str());
        }
        event.put(std::move(simedep));

        return;
    }
} // namespace
DEFINE_ART_MODULE(larg4::IonAndScint)
