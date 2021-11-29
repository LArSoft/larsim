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
// Restructured Nov. 21 by P. Green
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
#include "larsim/PhotonPropagation/PropagationTimeModel.h"

// Random number engine
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

#include <memory>
#include <vector>

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
  	bool 													fDoFastComponent;
    bool                          fDoSlowComponent;
    art::InputTag                 simTag;
    std::unique_ptr<ScintTime>    fScintTime;        // Tool to retrive timing of scintillation        
    CLHEP::HepRandomEngine&       fPhotonEngine;
    CLHEP::HepRandomEngine&       fScintTimeEngine;
    std::map<int, int>            PDChannelToSOCMapDirect; // Where each OpChan is.
    std::map<int, int>            PDChannelToSOCMapReflect; // Where each OpChan is.

    // propagation time model
    std::unique_ptr<PropagationTimeModel> fPropTimeModel;

    fhicl::ParameterSet fVUVTimingParams;
    fhicl::ParameterSet fVISTimingParams;

    bool fIncludePropTime;
  };
    
  //......................................................................    
  PDFastSimPVS::PDFastSimPVS(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fDoFastComponent(pset.get<bool>("DoFastComponent", true))
    , fDoSlowComponent(pset.get<bool>("DoSlowComponent", true))
    , simTag{pset.get<art::InputTag>("SimulationLabel")}
    , fScintTime{art::make_tool<ScintTime>(pset.get<fhicl::ParameterSet>("ScintTimeTool"))}
    , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "photon", pset, "SeedPhoton"))
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "scinttime", pset, "SeedScintTime"))
    , fIncludePropTime(pset.get<bool>("IncludePropTime", false))
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
    
    // Propagation times
    // validate configuration    
  	if(fIncludePropTime && !pset.get_if_present<fhicl::ParameterSet>("VUVTiming", fVUVTimingParams)) {
    throw art::Exception(art::errors::Configuration)
        << "Propagation time simulation requested, but VUVTiming not specified." << "\n";
  	}

  	if (pvs->StoreReflected() && fIncludePropTime && !pset.get_if_present<fhicl::ParameterSet>("VISTiming", fVISTimingParams)) {
      throw art::Exception(art::errors::Configuration)
          << "Reflected light propagation time simulation requested, but VISTiming not specified." << "\n";
    }

    // construct propagation time class
    if (fIncludePropTime) fPropTimeModel = std::make_unique<PropagationTimeModel>(fVUVTimingParams, fVISTimingParams, fScintTimeEngine, pvs->StoreReflected());

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
    for (unsigned int i = 0; i < nOpChannels; i ++) {    
			dir_photcol[i].fOpChannel  = i;
			ref_photcol[i].fOpChannel  = i;
			dir_phlitcol[i].OpChannel  = i;
			ref_phlitcol[i].OpChannel  = i;
		}
        
    art::Handle< std::vector<sim::SimEnergyDeposit> > edepHandle;
    if (!event.getByLabel(simTag, edepHandle)) {
			std::cout << "PDFastSimPVS Module Cannot getByLabel: " << simTag << std::endl;
			return;
    }
    else {
		std::cout << "PDFastSimPVS Module getByLabel: " << simTag << std::endl;
    }
        
    auto const& edeps = edepHandle;
        
    int num_points    = 0;
        
    for (auto const& edepi: *edeps) {
    	num_points ++;
            
			MappedCounts_t Visibilities;
			MappedCounts_t Visibilities_Ref;
            
			auto const& prt  = edepi.MidPoint();
			Visibilities     = pvs->GetAllVisibilities(prt);
			if(pvs->StoreReflected()) {
		    Visibilities_Ref = pvs->GetAllVisibilities(prt, true);
		    if(!Visibilities_Ref)	std::cout << "Fail to get visibilities for reflected photons." << std::endl;
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
			double edeposit   = edepi.Energy()/nphot;
			double pos[3]     = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
			int nphot_fast    = edepi.NumFPhotons();
			int nphot_slow    = edepi.NumSPhotons();

			// propagation time
		  std::vector<double> transport_time;

		  geo::Point_t const ScintPoint = {pos[0], pos[1], pos[2]};

			// loop through direct photons then reflected photons cases
		  for (size_t Reflected = 0; Reflected <= 1; ++Reflected) {

		  	// only do the reflected loop if including reflected light
		    if (Reflected && !pvs->StoreReflected()) continue;

		    // loop over each photo-detector
		    for (unsigned int channel = 0; channel < nOpChannels; ++ channel) {
		    	
		    	// visibility
		    	double visibleFraction;
		    	if (Reflected) visibleFraction = Visibilities_Ref[channel];   
		    	else visibleFraction = Visibilities[channel];

		    	if (visibleFraction == 0.0) continue; // voxel is not visible at this optical channel               

		    	// number of detected photons
		    	int ndetected_fast = 0;
			    int ndetected_slow = 0;

			    if (nphot_fast > 0) ndetected_fast = static_cast<int>(randpoisphot.fire(nphot_fast * visibleFraction));
			    if (nphot_slow > 0) ndetected_slow = static_cast<int>(randpoisphot.fire(nphot_slow * visibleFraction));

			    // calculate propagation times if included, does not matter whether fast or slow photon
		      if (fIncludePropTime) {
		      	transport_time.resize(ndetected_fast + ndetected_slow);
		        fPropTimeModel->propagationTime(transport_time, ScintPoint, channel, Reflected);
		      }

		      // SimPhotonsLite case
		      if (lgp->UseLitePhotons()) {

		      	sim::OpDetBacktrackerRecord tmpbtr(channel);

		      	if (ndetected_fast > 0 && fDoFastComponent) {
		      		for (long i = 0; i < ndetected_fast; ++i) {
								// calculate the time at which each photon is seen
								fScintTime->GenScintTime(true, fScintTimeEngine);
								int time;
								if (fIncludePropTime) time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + transport_time[i]);
								else time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
								if (Reflected) ++ref_phlitcol[channel].DetectedPhotons[time];
		            else ++dir_phlitcol[channel].DetectedPhotons[time];
								tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);                        
				      }
		      	}

		      	if ((ndetected_slow > 0) && fDoSlowComponent) {
				    	for (long i = 0; i < ndetected_slow; ++i) { 
				    		// calculate the time at which each photon is seen
								fScintTime->GenScintTime(false, fScintTimeEngine);
								int time;
		            if (fIncludePropTime) time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + transport_time[ndetected_fast + i]);
		            else time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
								if (Reflected) ++ref_phlitcol[channel].DetectedPhotons[time];
		            else ++dir_phlitcol[channel].DetectedPhotons[time];
								tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
							}
						}

						if (Reflected) AddOpDetBTR(*opbtr_ref, PDChannelToSOCMapReflect, tmpbtr);
						else AddOpDetBTR(*opbtr, PDChannelToSOCMapDirect, tmpbtr);		  
				  }

				  // SimPhotons case
				  else {

				  	sim::OnePhoton photon;
		        photon.SetInSD         = false;
		        photon.InitialPosition = edepi.End();
		        if (Reflected) photon.Energy = 2.9 * CLHEP::eV; // 430 nm
		        else photon.Energy = 9.7 * CLHEP::eV; // 128 nm

		        if (ndetected_fast > 0 && fDoFastComponent) {
		          for (long i = 0; i < ndetected_fast; ++i) {
		            // calculates the time at which the photon was produced
		            fScintTime->GenScintTime(true, fScintTimeEngine);
		            int time;
		            if (fIncludePropTime) time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + transport_time[i]);
		            else time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
		            photon.Time = time;
		            if(Reflected) ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
		            else dir_photcol[channel].insert(dir_photcol[channel].end(), 1, photon);
		          }
		        }

		        if (ndetected_slow > 0 && fDoSlowComponent) {
		          for (long i = 0; i < ndetected_slow; ++i) {
		            fScintTime->GenScintTime(false, fScintTimeEngine);
		            int time;
		            if (fIncludePropTime) time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + transport_time[ndetected_fast + i]);
		            else time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
		            photon.Time = time;
		            if(Reflected) ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
		            else dir_photcol[channel].insert(dir_photcol[channel].end(), 1, photon);
		          }
		        }
				  }
		    }
		  }
		}
  

	  PDChannelToSOCMapDirect.clear();
	  PDChannelToSOCMapReflect.clear();
	        
	  if (lgp->UseLitePhotons()) {
			event.put(move(phlit));
			event.put(move(opbtr));
			if (pvs->StoreReflected()) {
		    event.put(move(phlit_ref), "Reflected");
		    event.put(move(opbtr_ref), "Reflected");
		  }
	  }
	  else {
			event.put(move(phot));
			if (pvs->StoreReflected()) {
		    event.put(move(phot_ref), "Reflected");
		  }
	  }
	        
	  return;
	}
  /*                   
	for (unsigned int channel = 0; channel < nOpChannels; ++ channel)
	  {
	    auto visibleFraction     = Visibilities[channel];                
	    if (visibleFraction == 0.0)
	      {
		continue; //voxel is not visible at this optical channel.
	      }

	      // calculate number of photons visible at this channel
	                 
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
                                        
		if (pvs->StoreReflected() && Visibilities_Ref)
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
		photon.InitialPosition = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
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
                                        
		if (pvs->StoreReflected() && Visibilities_Ref)
		  {
		    auto visibleFraction_Ref = Visibilities_Ref [channel];
		    if (visibleFraction_Ref == 0.0)
		      {
			continue; //voxel is not visible at this optical channel.
		      }
		    sim::OnePhoton photon;
		    photon.SetInSD         = false;
		    photon.InitialPosition = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
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

      */
           
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
