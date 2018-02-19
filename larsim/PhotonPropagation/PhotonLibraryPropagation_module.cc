////////////////////////////////////////////////////////////////////////
// Class:       PhotonLibraryPropagation
// Plugin Type: producer (art v2_05_00)
// File:        PhotonLibraryPropagation_module.cc
//
// Generated at Tue Mar 21 07:45:42 2017 by Wesley Ketchum using cetskelgen
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
#include "nutools/RandomUtils/NuRandomService.h"

#include <memory>
#include <iostream>

#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/PhotonVoxels.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larsim/IonizationScintillation/ISCalculationSeparate.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


namespace phot {
  class PhotonLibraryPropagation;
}


class phot::PhotonLibraryPropagation : public art::EDProducer {
public:
  explicit PhotonLibraryPropagation(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonLibraryPropagation(PhotonLibraryPropagation const &) = delete;
  PhotonLibraryPropagation(PhotonLibraryPropagation &&) = delete;
  PhotonLibraryPropagation & operator = (PhotonLibraryPropagation const &) = delete;
  PhotonLibraryPropagation & operator = (PhotonLibraryPropagation &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob() override;

private:

  double fRiseTimeFast;
  double fRiseTimeSlow;
  bool   fDoSlowComponent;

  art::InputTag fEDepTag;

  larg4::ISCalculationSeparate fISAlg;

  double GetScintYield(sim::SimEnergyDeposit const&, detinfo::LArProperties const&);

  double GetScintTime(double scint_time, double rise_time, double, double);
  
  unsigned long N_VOXELS_X;
  unsigned long N_VOXELS_Y;
  unsigned long N_VOXELS_Z;
  std::unordered_map< unsigned long, std::vector< std::vector<double> > > fSCCalcMap;
  
  void PrecalculateSC(std::vector<sim::SimEnergyDeposit> const&);
  unsigned long GetSCMapIndex(float,float,float);
  std::vector<double> fPosOffsets;
};


phot::PhotonLibraryPropagation::PhotonLibraryPropagation(fhicl::ParameterSet const & p)
{
  produces< std::vector<sim::SimPhotons> >();
  art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "photon",    p, "SeedPhoton");
  art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "scinttime", p, "SeedScintTime");  
  this->reconfigure(p);
}

double phot::PhotonLibraryPropagation::GetScintYield(sim::SimEnergyDeposit const& edep,
						     detinfo::LArProperties const& larp)
{
  double yieldRatio = larp.ScintYieldRatio();
  if(larp.ScintByParticleType()){
    switch(edep.PdgCode()) {
    case 2212:
      yieldRatio = larp.ProtonScintYieldRatio();
      break;
    case 13:
    case -13:
      yieldRatio = larp.MuonScintYieldRatio();
      break;
    case 211:
    case -211:
      yieldRatio = larp.PionScintYieldRatio();
      break;
    case 321:
    case -321:
      yieldRatio = larp.KaonScintYieldRatio();
      break;
    case 1000020040:
      yieldRatio = larp.AlphaScintYieldRatio();
      break;
    case 11:
    case -11:
    case 22:
      yieldRatio = larp.ElectronScintYieldRatio();
      break;
    default:
      yieldRatio = larp.ElectronScintYieldRatio();
    }
  }
  return yieldRatio;
}

void phot::PhotonLibraryPropagation::produce(art::Event & e)
{
  art::ServiceHandle<PhotonVisibilityService> pvs;
  art::ServiceHandle<sim::LArG4Parameters> lgpHandle;
  const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();
  
  art::ServiceHandle<art::RandomNumberGenerator> rng;  
  CLHEP::HepRandomEngine &engine_photon = rng->getEngine("photon");
  CLHEP::RandPoissonQ randpoisphot(engine_photon);
  CLHEP::HepRandomEngine &engine_scinttime = rng->getEngine("scinttime");
  CLHEP::RandFlat randflatscinttime(engine_scinttime);
  
  const size_t NOpChannels = pvs->NOpChannels();
  double yieldRatio;
  int nphot,nphot_fast,nphot_slow;

  sim::OnePhoton photon;
  photon.Energy = 9.7e-6;
  photon.SetInSD = false;
  
  fISAlg.Initialize(larp,
		    lar::providerFrom<detinfo::DetectorPropertiesService>(),
		    &(*lgpHandle),
		    lar::providerFrom<spacecharge::SpaceChargeService>());
  
  auto const& edep_handle = e.getValidHandle< std::vector<sim::SimEnergyDeposit> >(fEDepTag);
  auto const& edeps(*edep_handle);

  std::unique_ptr< std::vector<sim::SimPhotons> > photCol ( new std::vector<sim::SimPhotons>);
  auto & photonCollection(*photCol);

  for(size_t i_op=0; i_op<NOpChannels; ++i_op){
    photonCollection.emplace_back(i_op);
    photonCollection[i_op].reserve(edeps.size()*10);
  }
  
  PrecalculateSC(edeps);
  
  for(auto const& edep : edeps){

    double const xyz[3] = { edep.X(), edep.Y(), edep.Z() };

    photon.InitialPosition = TVector3(xyz[0],xyz[1],xyz[2]);
    
    float const* Visibilities = pvs->GetAllVisibilities(xyz);
    if(!Visibilities)
      continue;
    
    yieldRatio = GetScintYield(edep,*larp);
    fISAlg.CalculateIonizationAndScintillation(edep,fSCCalcMap[GetSCMapIndex(edep.X(),edep.Y(),edep.Z())][1]);
    nphot =fISAlg.NumberScintillationPhotons();
    nphot_fast = yieldRatio*nphot;
    photon.Time = edep.T() + GetScintTime(larp->ScintFastTimeConst(),fRiseTimeFast,
					   randflatscinttime(),randflatscinttime());
    for(size_t i_op=0; i_op<NOpChannels; ++i_op)
      photonCollection[i_op].insert(photonCollection[i_op].end(),randpoisphot.fire(nphot_fast*Visibilities[i_op]),photon);

    if(fDoSlowComponent){
      nphot_slow = nphot - nphot_fast;
      
      if(nphot_slow>0){
	photon.Time = edep.T() - GetScintTime(larp->ScintSlowTimeConst(),fRiseTimeSlow,
					       randflatscinttime(),randflatscinttime());
	for(size_t i_op=0; i_op<NOpChannels; ++i_op)
	  photonCollection[i_op].insert(photonCollection[i_op].end(),randpoisphot.fire(nphot_slow*Visibilities[i_op]),photon);
      }
      
    }//end doing slow component

  }//end loop over edeps
  
  e.put(std::move(photCol));
  
}

double phot::PhotonLibraryPropagation::GetScintTime(double scint_time, double rise_time,
						    double r1, double r2)
{
  //no rise time
  if(rise_time<0.0)
    return -1 * scint_time * std::log(r1);

  while(1){
    double t = -1.0*scint_time*std::log(1-r1);
    double g = (scint_time+rise_time)/scint_time * std::exp(-1.0*t/scint_time)/scint_time;
    if ( r2 <= (std::exp(-1.0*t/rise_time)*(1-std::exp(-1.0*t/rise_time))/scint_time/scint_time*(scint_time+rise_time)) )
      return -1 * t;
  }
  
  return 1;
  
}

//----------------------------------------------------------------------------
unsigned long phot::PhotonLibraryPropagation::GetSCMapIndex(float x, float y, float z){
  unsigned long index = (unsigned long)x +
    N_VOXELS_X*((unsigned long)y) +
    N_VOXELS_X*N_VOXELS_Y*( (unsigned long)z );
  
  return index;
}
  
//----------------------------------------------------------------------------
void phot::PhotonLibraryPropagation::PrecalculateSC(std::vector<sim::SimEnergyDeposit> const& edeps)
{
  auto fSCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  art::ServiceHandle<geo::Geometry> geoHandle;
  auto const& geo = *geoHandle;

  fSCCalcMap.clear();
  for(auto const& edep : edeps){
    auto index = GetSCMapIndex(edep.X(),edep.Y(),edep.Z());
    if(fSCCalcMap[index].size()==0){
      double x = (double)((unsigned long)edep.X()) + 0.5;
      double y = (double)((unsigned long)edep.Y()) + 0.5 - (geo.DetHalfHeight());
      double z = (double)((unsigned long)edep.Z()) + 0.5;
      
	fSCCalcMap[index].push_back(fSCE->GetPosOffsets(edep.X(),edep.Y(),edep.Z()));
	fSCCalcMap[index].push_back(fSCE->GetEfieldOffsets(edep.X(),edep.Y(),edep.Z()));
    }
  }
  
}

void phot::PhotonLibraryPropagation::reconfigure(fhicl::ParameterSet const & p)
{
  fRiseTimeFast = p.get<double>("RiseTimeFast",-1.0);
  fRiseTimeSlow = p.get<double>("RiseTimeSlow",-1.0);
  fDoSlowComponent = p.get<bool>("DoSlowComponent");

  fEDepTag = p.get<art::InputTag>("EDepModuleLabel");
}

void phot::PhotonLibraryPropagation::beginJob()
{

  art::ServiceHandle<geo::Geometry> geoHandle;
  auto const& geo = *geoHandle;

  N_VOXELS_X = (unsigned long)(geo.DetHalfWidth()*2) + 1;
  N_VOXELS_Y = (unsigned long)(geo.DetHalfHeight()*2) + 1;
  N_VOXELS_Z = (unsigned long)(geo.DetLength()) + 1;

  
}

DEFINE_ART_MODULE(phot::PhotonLibraryPropagation)
