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
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"
#include "larsim/Simulation/LArG4Parameters.h"

// Random number engine
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
//#include "CLHEP/Random/RandGauss.h"

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
        bool                          fDoSlowComponent;
        art::InputTag                 simTag;
        std::unique_ptr<ScintTime>    fScintTime;        // Tool to retrive timinig of scintillation        
        CLHEP::HepRandomEngine&       fPhotonEngine;
        CLHEP::HepRandomEngine&       fScintTimeEngine;
        std::map<int, int>            PDChannelToSOCMapDirect; // Where each OpChan is.
        std::map<int, int>            PDChannelToSOCMapReflect; // Where each OpChan is.
    };
    
    //......................................................................    
    PDFastSimPVS::PDFastSimPVS(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fDoSlowComponent{pset.get<bool>("DoSlowComponent")}
    , simTag{pset.get<art::InputTag>("SimulationLabel")}
    , fScintTime{art::make_tool<ScintTime>(pset.get<fhicl::ParameterSet>("ScintTimeTool"))}
    , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "photon", pset, "SeedPhoton"))
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "scinttime", pset, "SeedScintTime"))
    {
        std::cout << "PDFastSimPVS Module Construct" << std::endl;
        
        art::ServiceHandle<PhotonVisibilityService const> pvs;
        art::ServiceHandle<sim::LArG4Parameters const> lgp;        
        if (lgp->UseLitePhotons())
        {
            std::cout << "Use Lite Photon." << std::endl;
            produces< std::vector<sim::SimPhotonsLite> >();
            produces< std::vector<sim::OpDetBacktrackerRecord> >();
            
            if(pvs->StoreReflected())
            {
                std::cout << "Store Reflected Photons" << std::endl;
                produces< std::vector<sim::SimPhotonsLite> >("Reflected");
                produces< std::vector<sim::OpDetBacktrackerRecord> >("Reflected");     
            }            
        }
        else
        {
            std::cout << "Use Sim Photon." << std::endl;
            produces< std::vector<sim::SimPhotons> >();
            if(pvs->StoreReflected())
            {
                std::cout << "Store Reflected Photons" << std::endl;            
                produces< std::vector<sim::SimPhotons> >("Reflected");     
            }            
        }        
    }
    
    //......................................................................    
    void PDFastSimPVS::produce(art::Event& event)
    {
        std::cout << "PDFastSimPVS Module Producer" << std::endl;
        
        art::ServiceHandle<PhotonVisibilityService const> pvs;
        art::ServiceHandle<sim::LArG4Parameters const> lgp;
        auto const nOpChannels = pvs->NOpChannels();
        
        CLHEP::RandPoissonQ randpoisphot{fPhotonEngine};
        
        std::unique_ptr< std::vector< sim::SimPhotons > >             phot   (new std::vector<sim::SimPhotons>);
        std::unique_ptr< std::vector< sim::SimPhotonsLite > >         phlit  (new std::vector<sim::SimPhotonsLite>);
        std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > > opbtr  (new std::vector<sim::OpDetBacktrackerRecord>);
        
        std::unique_ptr< std::vector< sim::SimPhotons > >             phot_ref   (new std::vector<sim::SimPhotons>);
        std::unique_ptr< std::vector< sim::SimPhotonsLite > >         phlit_ref  (new std::vector<sim::SimPhotonsLite>);
        std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > > opbtr_ref  (new std::vector<sim::OpDetBacktrackerRecord>);
        
        auto& dir_photcol(*phot);
        auto& ref_photcol(*phot_ref);
        auto& dir_phlitcol(*phlit);
        auto& ref_phlitcol(*phlit_ref);
        dir_photcol.resize(nOpChannels);
        ref_photcol.resize(nOpChannels);
        dir_phlitcol.resize(nOpChannels);
        ref_phlitcol.resize(nOpChannels);
        for (unsigned int i = 0; i < nOpChannels; i ++)
        {
            dir_photcol[i].fOpChannel  = i;
            ref_photcol[i].fOpChannel  = i;
            dir_phlitcol[i].OpChannel  = i;
            ref_phlitcol[i].OpChannel  = i;
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
        
        for (auto const& edepi: *edeps)
        {
            num_points ++;
            
            MappedCounts_t Visibilities;
            MappedCounts_t Visibilities_Ref;
            
            auto const& prt  = edepi.MidPoint();
            Visibilities     = pvs->GetAllVisibilities(prt);
            if(pvs->StoreReflected())
            {
                Visibilities_Ref = pvs->GetAllVisibilities(prt, true);
            }
            
            if(!Visibilities)
            {
                //throw cet::exception("PDFastSimPVS")
                std::cout << "There is no entry in the PhotonLibrary for this position in space. Position: " << edepi.MidPoint();
                std::cout << "\n Move to next point" << std::endl;
                continue;
            }
            int trackID       = edepi.TrackID();
            int nphot         = edepi.NumPhotons();
            int nion          = edepi.NumElectrons();
            double edeposit   = edepi.Energy()/nphot;
            double pos[3]     = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
            
            int nphot_fast    = edepi.NumFPhotons();
            int nphot_slow    = edepi.NumSPhotons();
                        
            for (unsigned int channel = 0; channel < nOpChannels; ++ channel)
            {
                auto visibleFraction     = Visibilities[channel];                
                if (visibleFraction == 0.0)
                {
                    continue; //voxel is not visible at this optical channel.
                }
                
                if (lgp->UseLitePhotons())
                {
                    sim::OpDetBacktrackerRecord tmpbtr(channel);
                    if (nphot_fast > 0)
                    {
                        //random number, poisson distribution, mean: the amount of photons visible at this channel
                        auto n = static_cast<int>(randpoisphot.fire(nphot_fast * visibleFraction));
                        for (long i = 0; i < n; ++i) 
                        {
                            //calculates the time at which the photon was produced
                            fScintTime->GenScintTime(true, fScintTimeEngine);
                            auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                            ++ dir_phlitcol[channel].DetectedPhotons[time];
                            tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);                        
                        }
                    }
                    
                    if ((nphot_slow > 0) && fDoSlowComponent) 
                    {
                        auto n = static_cast<int>(randpoisphot.fire(nphot_slow * visibleFraction));
                        for (long i = 0; i < n; ++i) 
                        {
                            fScintTime->GenScintTime(false, fScintTimeEngine);
                            auto time = static_cast<int>(edepi.StartT()+ fScintTime->GetScintTime());
                            ++ dir_phlitcol[channel].DetectedPhotons[time];
                            tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);                    }
                    }
                    
                    AddOpDetBTR(*opbtr, PDChannelToSOCMapDirect, tmpbtr);

                    if (pvs->StoreReflected())
                    {
                        sim::OpDetBacktrackerRecord tmpbtr_ref(channel); 
                        auto visibleFraction_Ref = Visibilities_Ref [channel];                        
                        if (visibleFraction_Ref == 0.0)
                        {
                            continue; //voxel is not visible at this optical channel.
                        }
                        if (nphot_fast > 0)
                        {
                            auto n = static_cast<int>(randpoisphot.fire(nphot_fast * visibleFraction_Ref));
                            for (long i = 0; i < n; ++i)
                            {
                                //calculates the time at which the photon was produced
                                fScintTime->GenScintTime(true, fScintTimeEngine);
                                auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                                ++ ref_phlitcol[channel].DetectedPhotons[time];
                                tmpbtr_ref.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
                            }
                        }
                        
                        if ((nphot_slow > 0) && fDoSlowComponent)
                        {
                            auto n = static_cast<int>(randpoisphot.fire(nphot_slow * visibleFraction_Ref));
                            for (long i = 0; i < n; ++i)
                            {
                                fScintTime->GenScintTime(false, fScintTimeEngine);
                                auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                                ++ ref_phlitcol[channel].DetectedPhotons[time];
                                tmpbtr_ref.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
                            }
                        }
                        
                        AddOpDetBTR(*opbtr_ref, PDChannelToSOCMapReflect, tmpbtr_ref);
                    }
                }
                else
                {
                    sim::OnePhoton photon;
                    photon.SetInSD         = false;
                    photon.InitialPosition = edepi.End();
                    photon.Energy          = 9.7e-6;

                    if (nphot_fast > 0)
                    {
                        //random number, poisson distribution, mean: the amount of photons visible at this channel
                        auto n = static_cast<int>(randpoisphot.fire(nphot_fast * visibleFraction));
                        if (n > 0)
                        {
                            fScintTime->GenScintTime(true, fScintTimeEngine);
                            auto time   = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                            photon.Time = time;
                            // add n copies of sim::OnePhoton photon to the photon collection for a given OpChannel
                            dir_photcol[channel].insert(dir_photcol[channel].end(), n, photon);
                        }
                    }
                    if ((nphot_slow > 0) && fDoSlowComponent)
                    {
                        //random number, poisson distribution, mean: the amount of photons visible at this channel
                        auto n = static_cast<int>(randpoisphot.fire(nphot_slow * visibleFraction));
                        if (n > 0)
                        {
                            fScintTime->GenScintTime(false, fScintTimeEngine);
                            auto time   = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                            photon.Time = time;
                            // add n copies of sim::OnePhoton photon to the photon collection for a given OpChannel
                            dir_photcol[channel].insert(dir_photcol[channel].end(), n, photon);
                        }
                    }
                    
                    if(pvs->StoreReflected())
                    {
                        auto visibleFraction_Ref = Visibilities_Ref [channel];
                        if (visibleFraction_Ref == 0.0)
                        {
                            continue; //voxel is not visible at this optical channel.
                        }
                        sim::OnePhoton photon;
                        photon.SetInSD         = false;
                        photon.InitialPosition = edepi.End();
                        photon.Energy          = 2.7e-6;
                        
                        if (nphot_fast > 0)
                        {
                            auto n = static_cast<int>(randpoisphot.fire(nphot_fast * visibleFraction_Ref));
                            for (long i = 0; i < n; ++i)
                            {
                                fScintTime->GenScintTime(true, fScintTimeEngine);
                                auto time   = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                                photon.Time = time;
                                ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
                            }
                        }
                        
                        if ((nphot_slow > 0) && fDoSlowComponent)
                        {
                            auto n = static_cast<int>(randpoisphot.fire(nphot_slow * visibleFraction_Ref));
                            for (long i = 0; i < n; ++i)
                            {
                                fScintTime->GenScintTime(false, fScintTimeEngine);
                                auto time   = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                                photon.Time = time;
                                ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
                            }
                        }
                    }
                }
            }
        }
        
        PDChannelToSOCMapDirect.clear();
        PDChannelToSOCMapReflect.clear();
        
        if (lgp->UseLitePhotons())
        {
            event.put(move(phlit));
            event.put(move(opbtr));
            if (pvs->StoreReflected())
            {
                event.put(move(phlit_ref), "Reflected");
                event.put(move(opbtr_ref), "Reflected");
            }
        }
        else
        {
            event.put(move(phot));
            if (pvs->StoreReflected())
            {
                event.put(move(phot_ref), "Reflected");
            }
        }
        
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
