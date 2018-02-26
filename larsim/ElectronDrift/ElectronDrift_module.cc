////////////////////////////////////////////////////////////////////////
// Class:       ElectronDrift
// Plugin Type: producer (art v2_05_00)
// File:        ElectronDrift_module.cc
//
// Generated at Sat Mar 11 17:31:14 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include <memory>
#include <unordered_map>
#include <chrono>

#include "nutools/RandomUtils/NuRandomService.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcore/Geometry/Geometry.h"

#include "larsim/ElectronDrift/ElectronDriftAlg.hh"

namespace larg4 {
  class ElectronDrift;
}


class larg4::ElectronDrift : public art::EDProducer {
public:
  explicit ElectronDrift(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ElectronDrift(ElectronDrift const &) = delete;
  ElectronDrift(ElectronDrift &&) = delete;
  ElectronDrift & operator = (ElectronDrift const &) = delete;
  ElectronDrift & operator = (ElectronDrift &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob() override;

private:

  ElectronDriftAlg       fElectronDriftAlg;
  art::InputTag          fEDepTag;

  //map to get position in simChannelCol of channelID
  std::unordered_map<raw::ChannelID_t,size_t> fChannelMap;

  std::vector<double>    fRecipDriftVelocities;
  
};


larg4::ElectronDrift::ElectronDrift(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
  art::ServiceHandle<rndm::NuRandomService>()
    ->createEngine(*this, "HepJamesRandom", "ElectronDrift", p, "Seed");

  produces< std::vector<sim::SimChannel> >();

}

void larg4::ElectronDrift::produce(art::Event & e)
{

  std::unique_ptr< std::vector<sim::SimChannel> > simChannelCol(new std::vector<sim::SimChannel>);
  auto & simch_vec(*simChannelCol);

  //some variables we will use (could be private, but probably doesn't matter)
  double XYZ[3];
  double edep_XYZ[3];
  geo::TPCID tpc_id;
  unsigned int tdc;
  raw::ChannelID_t channel;
  double tdrift;
  std::unordered_map<raw::ChannelID_t,size_t>::iterator it_chmap;
  
  //services
  art::ServiceHandle<geo::Geometry> geoHandle;
  auto const& geo = *geoHandle;
  auto const& detclk = (*lar::providerFrom<detinfo::DetectorClocksService>());
  auto const& detprop = (*lar::providerFrom<detinfo::DetectorPropertiesService>());

  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine("ElectronDrift");
  
    //some setup using the services...
  art::ServiceHandle<sim::LArG4Parameters> lgpHandle;
  fElectronDriftAlg.Initialize(lar::providerFrom<detinfo::LArPropertiesService>(),
			       &detprop,
			       &(*lgpHandle),
			       lar::providerFrom<spacecharge::SpaceChargeService>(),
			       &geo,
			       &engine);
  
  
  fRecipDriftVelocities.resize(geo.Nplanes());
  for (size_t i = 0; i<geo.Nplanes(); ++i)
    fRecipDriftVelocities[i] = 1./(detprop.DriftVelocity(detprop.Efield(i),
							 detprop.Temperature())/1000.);

  
  //ok, grab the edeps
  auto const& edep_handle = e.getValidHandle< std::vector<sim::SimEnergyDeposit> >(fEDepTag);
  auto const& edep_vec(*edep_handle);

  fChannelMap.clear();
  simch_vec.reserve(geo.Nchannels());

  double time_drift = 0;
  double time_simch = 0;

  //auto t00 = std::chrono::high_resolution_clock::now();
  fElectronDriftAlg.PrecalculateSC(edep_vec);
  //auto t0 = std::chrono::high_resolution_clock::now();
  
  //std::chrono::duration<double, std::milli> tmp0;
  //tmp0 = t0 - t00;
  //std::cout << "Precalc time was " << tmp0.count() << std::endl;

  //loop over edeps
  for(auto const& edep : edep_vec){

    //run the drifiting...
    //auto t1 = std::chrono::high_resolution_clock::now();
    fElectronDriftAlg.DriftElectrons(edep);
    //auto t2 = std::chrono::high_resolution_clock::now();
    
    edep_XYZ[0] = edep.X(); edep_XYZ[1] = edep.Y(); edep_XYZ[2] = edep.Z();
    
    XYZ[0] = edep.X(); //we won't need to validate x position: that will show up as timing diff, that's all

    /*
    std::cout << "\tDrifted electrons for edep at " 
	      << edep.X() << "," << edep.Y() << "," << edep.Z()
	      << " is " << fElectronDriftAlg.NElectronClusters() << std::endl;
    */
    //loop over the electron clusters and add them to channels
    for(size_t i_cl=0; i_cl < fElectronDriftAlg.NElectronClusters(); ++i_cl){
      
      //check to make sure diffused cluster is inside the original tpc
      XYZ[1] = fElectronDriftAlg.ElectronClustersY()[i_cl];
      XYZ[2] = fElectronDriftAlg.ElectronClustersZ()[i_cl];

      tpc_id = geo.FindTPCAtPosition(XYZ);
      if(tpc_id.TPC==geo::TPCID::InvalidID)
	continue;
      auto const& tpcgeo = geo.TPC(tpc_id);

      tdrift = fElectronDriftAlg.ElectronClustersT()[i_cl];
      
      //ok, so, per plane, let's add the clusters to each channel
      for(size_t i_p=0; i_p < tpcgeo.Nplanes(); ++i_p){

	//increment the drift time for the next plane if this is not first plane
	if(i_p!=0) tdrift += tpcgeo.PlanePitch(i_p-1,i_p) * fRecipDriftVelocities[i_p];
		     
	XYZ[0] = tpcgeo.PlaneLocation(0)[0] - tpcgeo.Plane0Pitch(i_p);
	try{
	channel = geo.NearestChannel(XYZ,i_p,tpc_id.TPC,tpc_id.Cryostat);
	tdc     = detclk.TPCClock().Ticks(detclk.G4ToElecTime(tdrift));
	
	it_chmap = fChannelMap.find(channel);
	size_t loc=0;
	if(it_chmap==fChannelMap.end()){
	  loc = simch_vec.size();
	  simch_vec.emplace_back(channel);
	  fChannelMap[channel] = simch_vec.size()-1;
	}
	else{
	  loc = it_chmap->second;
	}
	simch_vec[loc].AddIonizationElectrons(edep.TrackID(),
					      tdc,
					      fElectronDriftAlg.ElectronClustersEl()[i_cl],
					      edep_XYZ,
					      fElectronDriftAlg.ElectronClustersEn()[i_cl]);
	}
	catch(...){}
      }//end loop over planes
    }//end loop over electron clusters

    //auto t3 = std::chrono::high_resolution_clock::now();

    //std::chrono::duration<double, std::milli> tmp;
    //tmp = t2 - t1; time_drift += tmp.count();
    //tmp = t3 - t2; time_simch += tmp.count();
    
  }//end loop over energy depositions

  //std::cout << "TOTAL ELAPSED TIME IS " << time_drift << " FOR DRIFT AND " << time_simch << " FOR SIMCH" << std::endl;

  //for(auto const& time : fElectronDriftAlg.GetTimer())
  //std::cout << "\tEdep time check ... " << time << std::endl;
  
  //put sim channles onto the event
  std::cout << "Going to put " << simChannelCol->size() << " channels on the event." << std::endl;
  e.put(std::move(simChannelCol));
}

void larg4::ElectronDrift::reconfigure(fhicl::ParameterSet const & p)
{
  fEDepTag = p.get<art::InputTag>("EDepModuleLabel");
}

void larg4::ElectronDrift::beginJob()
{
}

DEFINE_ART_MODULE(larg4::ElectronDrift)
