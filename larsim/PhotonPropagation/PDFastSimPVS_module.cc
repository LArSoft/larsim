////////////////////////////////////////////////////////////////////////
// Class:       PDFastSimPVS
// Plugin Type: producer
// File:        PDFastSimPVS_module.cc
// Description:
// - acts on sim::SimEnergyDeposit from LArG4Main, 
// - simulate (fast, photon visibility service) the OpDet response to optical photons
// Input: 'sim::SimEnergyDeposit'
// Output: 'sim::OpDetBacktrackerRecord'
//Fast simulation of propagating the photons created from SimEnergyDeposits.

//This module does a fast simulation of propagating the photons created from SimEnergyDeposits,
//This simulation is done using the PhotonLibrary, which stores the visibilities of each optical channel
//with respect to each optical voxel in the TPC volume, to avoid propagating single photons using Geant4.
//At the end of this module a collection of the propagated photons either as
//'sim::OpDetBacktrackerRecord' are placed into the art event.

//The steps this module takes are:
//  - to take number of photon and the vertex information from 'sim::SimEnergyDeposits',
//  - use the PhotonLibrary (visibilities) to determine the amount of visible photons at each optical channel,
//  - visible photons: the number of photons times the visibility at the middle of the Geant4 step for a given optical channel.
//  - other photon information is got from 'sim::SimEnergyDeposits'
//  - add 'sim::OpDetBacktrackerRecord' to event
// Aug. 19 by Mu Wei
////////////////////////////////////////////////////////////////////////

// Art libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

// Random number engine
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
//#include "CLHEP/Random/RandGauss.h"

using namespace std;

namespace
{
    
    //......................................................................    
    double single_exp(double t, double tau2)
    {
        return exp((-1.0 * t) / tau2) / tau2;
    }

    double bi_exp(double t, double tau1, double tau2)
    {
        return (((exp((-1.0 * t) / tau2) * (1.0 - exp((-1.0 * t) / tau1))) / tau2) / tau2) * (tau1 + tau2);
    }

    
    //......................................................................    
    // Returns the time within the time distribution of the scintillation process, when the photon was created.
    // Scintillation light has an exponential decay which is given by the decay time, tau2,
    // and an exponential increase, which here is given by the rise time, tau1.
    // randflatscinttime is passed to use the saved seed from the RandomNumberSaver in order to be able to reproduce the same results.
    double GetScintTime(double tau1, double tau2, CLHEP::RandFlat& randflatscinttime)
    {
        // tau1: rise time (originally defaulted to -1) and tau2: decay time
        //ran1, ran2 = random numbers for the algorithm
        if ((tau1 == 0.0) || (tau1 == -1.0))
        {
            return -tau2 * log(randflatscinttime());
        }
        while (1)
        {
            auto ran1 = randflatscinttime();
            auto ran2 = randflatscinttime();
            auto d = (tau1 + tau2) / tau2;
            auto t = -tau2 * log(1 - ran1);
            auto g = d * single_exp(t, tau2);
            if (ran2 <= bi_exp(t, tau1, tau2) / g)
            {
                return t;
            }
        }
    }
} // namespace std

namespace phot
{
    class PDFastSimPVS : public art::EDProducer
    {
    public:
        explicit PDFastSimPVS(fhicl::ParameterSet const&);
        void produce(art::Event&) override;
        void AddOpDetBTR(std::vector< sim::OpDetBacktrackerRecord > & opbtr,
                             std::map<int, int> & ChannelMap,
                             sim::OpDetBacktrackerRecord btr);
    private:
        double                  fRiseTimeFast;
        double                  fRiseTimeSlow;
        bool                    fDoSlowComponent;
        art::InputTag           simTag;
        CLHEP::HepRandomEngine& fPhotonEngine;
        CLHEP::HepRandomEngine& fScintTimeEngine;
        std::map<int, int>      PDChannelToSOCMap; //Where each OpChan is.
    };
    
    //......................................................................    
    PDFastSimPVS::PDFastSimPVS(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fRiseTimeFast{pset.get<double>("RiseTimeFast", 0.0)}
    , fRiseTimeSlow{pset.get<double>("RiseTimeSlow", 0.0)}
    , fDoSlowComponent{pset.get<bool>("DoSlowComponent")}
    , simTag{pset.get<art::InputTag>("SimulationLabel")}
    , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "photon", pset, "SeedPhoton"))
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "scinttime", pset, "SeedScintTime"))
    {
        std::cout << "PDFastSimPVS Module Construct" << std::endl;
        
        produces< std::vector<sim::SimPhotonsLite> >("pvs");
        produces< std::vector<sim::OpDetBacktrackerRecord> >("pvs");     
    }
    
    //......................................................................    
    void PDFastSimPVS::produce(art::Event& event)
    {
        std::cout << "PDFastSimPVS Module Producer" << std::endl;
        
        art::ServiceHandle<PhotonVisibilityService const> pvs;
        auto const* larp = lar::providerFrom<detinfo::LArPropertiesService>();
        auto const nOpChannels = pvs->NOpChannels();
        
        CLHEP::RandPoissonQ randpoisphot{fPhotonEngine};
        CLHEP::RandFlat randflatscinttime{fScintTimeEngine};
        
        std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > > opbtr  (new std::vector<sim::OpDetBacktrackerRecord>);
        std::unique_ptr< std::vector< sim::SimPhotonsLite> >          phlit  (new std::vector<sim::SimPhotonsLite>);
        
        auto& photonLiteCollection (*phlit);
        photonLiteCollection.resize(nOpChannels);
        for (unsigned int i = 0; i < nOpChannels; i ++)
        {
            photonLiteCollection[i].OpChannel = i;
        }
        
        art::Handle< std::vector<sim::SimEnergyDeposit> > edepHandle;
        if (!event.getByLabel(simTag, edepHandle))
        {
            std::cout << "PDFastSimPVS Module Cannot getByLabel: " << simTag << std::endl;
            return;
        }
                
        art::ServiceHandle<geo::Geometry> geom;
        auto const& edeps = edepHandle;
        
        int num_points    = 0;
        int num_fastph    = 0;
        int num_slowph    = 0;
        int num_fastdp    = 0;
        int num_slowdp    = 0;
        
        float vis_scale   = 1.0; // to scale the visibility fraction, for test only;
        
        for (auto const& edepi: *edeps)
        {
            num_points ++;
            
            auto const& prt          = edepi.MidPoint();
            auto const& Visibilities = pvs->GetAllVisibilities(prt);
            if (!Visibilities)
            {
                //throw cet::exception("PDFastSimPVS")
                std::cout << "There is no entry in the PhotonLibrary for this position in space. "
                "Position: " << edepi.MidPoint();
                std::cout << "\n Move to next point" << std::endl;
                continue;
            }
            
            int trackID       = edepi.TrackID();
            double nphot      = edepi.NumPhotons();
            double edeposit   = edepi.Energy()/nphot;
            double pos[3]     = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
            
            double nphot_fast = edepi.NumFPhotons();
            double nphot_slow = edepi.NumSPhotons();
            
            num_fastph += nphot_fast;
            num_slowph += nphot_slow;
            
            for (unsigned int channel = 0; channel < nOpChannels; ++ channel)
            {
                auto visibleFraction = Visibilities[channel] * vis_scale;  // to scale the visibility fraction;
                
                if (visibleFraction == 0.0)
                {
                    continue; //voxel is not visible at this optical channel.
                }
                
                sim::OpDetBacktrackerRecord tmpbtr(channel);
                
                if (nphot_fast > 0)
                {
                    //random number, poisson distribution, mean: the amount of photons visible at this channel
                    auto n = static_cast<int>(randpoisphot.fire(nphot_fast * visibleFraction));
                    num_fastdp += n;                    
                    for (long i = 0; i < n; ++i) 
                    {
                        //calculates the time at which the photon was produced
                        auto time = static_cast<int>(edepi.StartT() + GetScintTime(fRiseTimeFast, larp->ScintFastTimeConst(), randflatscinttime));
                        ++ photonLiteCollection[channel].DetectedPhotons[time];
                        tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);                        
                    }
                }
                
                if ((nphot_slow > 0) && fDoSlowComponent) 
                {
                    auto n = static_cast<int>(randpoisphot.fire(nphot_slow * visibleFraction));
                    num_slowdp += n;                    
                    for (long i = 0; i < n; ++i) 
                    {
                        auto time = static_cast<int>(edepi.StartT() + GetScintTime(fRiseTimeSlow, larp->ScintSlowTimeConst(), randflatscinttime));
                        ++ photonLiteCollection[channel].DetectedPhotons[time];
                        tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);                    }
                }
                
                AddOpDetBTR(*opbtr, PDChannelToSOCMap, tmpbtr);
            }
        }
        
        std::cout << "Total points: " << num_points << ", total fast photons: " << num_fastph << ", total slow photons: " << num_slowph << std::endl;
        std::cout << "detected fast photons: " << num_fastdp << ", detected slow photons: " << num_slowdp << std::endl;
        
        PDChannelToSOCMap.clear();
        event.put(move(phlit), "pvs");
        event.put(move(opbtr), "pvs");
        
        return;
    }
    
    //......................................................................    
    void PDFastSimPVS::AddOpDetBTR(std::vector< sim::OpDetBacktrackerRecord > & opbtr,
                                         std::map<int, int> & ChannelMap,
                                         sim::OpDetBacktrackerRecord btr) 
    {
        int iChan = btr.OpDetNum();
        std::map<int, int>::iterator channelPosition = ChannelMap.find(iChan);
        
        if (channelPosition == ChannelMap.end() )
        {
            ChannelMap[iChan] = opbtr.size();
            opbtr.emplace_back(std::move(btr));
        }
        else
        {
            unsigned int idtest = channelPosition->second;
            auto const& timePDclockSDPsMap = btr.timePDclockSDPsMap();
            
            for(auto const& timePDclockSDP : timePDclockSDPsMap)
            {
                for(auto const& sdp : timePDclockSDP.second)
                {
                    double xyz[3] = {sdp.x, sdp.y, sdp.z};
                    opbtr.at(idtest).AddScintillationPhotons(sdp.trackID,
                                                             timePDclockSDP.first,
                                                             sdp.numPhotons,
                                                             xyz,
                                                             sdp.energy);
                }
            }
        }    
    }  
} // namespace

DEFINE_ART_MODULE(phot::PDFastSimPVS)
