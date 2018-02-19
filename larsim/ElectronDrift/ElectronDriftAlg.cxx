////////////////////////////////////////////////////////////////////////
/// \file  ElectronDrift.hh
/// \brief Algorithm class for modeling electron drift, given an
///        energy deposition in space.
///
/// \version $Id:  $
/// \author  wketchum@fnal.gov
////////////////////////////////////////////////////////////////////////
#include "larsim/ElectronDrift/ElectronDriftAlg.hh"

#include "lardata/DetectorInfo/LArProperties.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larcore/Geometry/GeometryCore.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "TMath.h"

#include <random>

namespace larg4{

  //----------------------------------------------------------------------------
  ElectronDriftAlg::ElectronDriftAlg()
  {
  }

  //----------------------------------------------------------------------------
  ElectronDriftAlg::~ElectronDriftAlg()
  {
  }

  //----------------------------------------------------------------------------
  void ElectronDriftAlg::Initialize(const detinfo::LArProperties* larp,   
				    const detinfo::DetectorProperties* detp,
				    const sim::LArG4Parameters* lgp,
				    const spacecharge::SpaceCharge* sce,
				    const geo::GeometryCore* geo,
				    CLHEP::HepRandomEngine* engine)
  {
    fSCE = sce;
    fGeo = geo;
    fEngine = engine;
    
    fElectronLifetime      = -1000. * detp->ElectronLifetime();
    fDriftVelocityVolume   = detp->DriftVelocity(detp->Efield(0),
						 detp->Temperature())/1000.;
    fRecipDriftVelocityVolume = 1. / fDriftVelocityVolume;
    fLongitudinalDiffConst = lgp->LongitudinalDiffusion();
    fTransverseDiffConst   = lgp->TransverseDiffusion();

    fElectronClusterSize   = lgp->ElectronClusterSize();
    fMinNumberOfElCluster  = lgp->MinNumberOfElCluster();

    fISAlg.Initialize(larp,detp,lgp,sce);

    N_VOXELS_X = (unsigned long)(fGeo->DetHalfWidth()*2) + 1;
    N_VOXELS_Y = (unsigned long)(fGeo->DetHalfHeight()*2) + 1;
    N_VOXELS_Z = (unsigned long)(fGeo->DetLength()) + 1;

    timer.clear();
    timer.resize(8,0.0);
  }

  //----------------------------------------------------------------------------
  unsigned long ElectronDriftAlg::GetSCMapIndex(float x, float y, float z){
    unsigned long index = (unsigned long)x +
      N_VOXELS_X*((unsigned long)y) +
      N_VOXELS_X*N_VOXELS_Y*( (unsigned long)z );

    return index;
  }
  
  //----------------------------------------------------------------------------
  void ElectronDriftAlg::PrecalculateSC(std::vector<sim::SimEnergyDeposit> const& edeps)
  {
    fSCCalcMap.clear();
    for(auto const& edep : edeps){
      auto index = GetSCMapIndex(edep.X(),edep.Y(),edep.Z());
      if(fSCCalcMap[index].size()==0){
	double x = (double)((unsigned long)edep.X()) + 0.5;
	double y = (double)((unsigned long)edep.Y()) + 0.5 - (fGeo->DetHalfHeight());
	double z = (double)((unsigned long)edep.Z()) + 0.5;

	fSCE->GetPosOffsets(edep.X(),edep.Y(),edep.Z(),fPosOffsets);	
	fSCCalcMap[index].push_back(fPosOffsets);
	fSCE->GetEfieldOffsets(edep.X(),edep.Y(),edep.Z(),fPosOffsets);	
	fSCCalcMap[index].push_back(fPosOffsets);
      }
    }

  }
  

  //----------------------------------------------------------------------------
  void ElectronDriftAlg::DriftElectrons(sim::SimEnergyDeposit const& edep)
  {
    //auto t1 = std::chrono::high_resolution_clock::now();
    fNElectronClusters = 0;
    CLHEP::RandGaussQ randGauss(*fEngine);
    std::random_device rd;
    std::mt19937 gen(rd());
    //auto t2 = std::chrono::high_resolution_clock::now();

    //get the right TPC
    fXYZ[0] = edep.X(); fXYZ[1] = edep.Y(); fXYZ[2] = edep.Z(); 
    fTPC_id = fGeo->FindTPCAtPosition(fXYZ);
    if(fTPC_id.TPC==geo::TPCID::InvalidID){
      //throw cet::exception("DriftElectrons")
      std::cout << "Invalid TPC ID for position "
		<< edep.X() << " , " << edep.Y() << " , " << edep.Z()
		<< " --- " 
		<< fXYZ[0] << " , " << fXYZ[1] << " , " << fXYZ[2]
		<< std::endl;
      return;
    }
    auto const& tpcgeo = fGeo->TPC(fTPC_id);
    //auto t3 = std::chrono::high_resolution_clock::now();

    fPosOffsets = fSCCalcMap[GetSCMapIndex(edep.X(),edep.Y(),edep.Z())][0];
    //auto t4 = std::chrono::high_resolution_clock::now();

    //determine xdrift distance
    if(tpcgeo.DriftDirection() == geo::kNegX)
      fXDrift = edep.X() - tpcgeo.PlaneLocation(0)[0] - fPosOffsets[0];
    else if(tpcgeo.DriftDirection() == geo::kNegX)
      fXDrift = tpcgeo.PlaneLocation(0)[0] - edep.X() - fPosOffsets[0];
    if(fXDrift<0)
      return;

    //get drift time
    fTDrift = fXDrift * fRecipDriftVelocityVolume;
    //auto t5 = std::chrono::high_resolution_clock::now();

    //get n_electrons
    fISAlg.Reset();
    fISAlg.CalculateIonizationAndScintillation(edep,fSCCalcMap[GetSCMapIndex(edep.X(),edep.Y(),edep.Z())][1]);
    fNElectrons = fISAlg.NumberIonizationElectrons() * TMath::Exp(fTDrift / fElectronLifetime);
    if(fNElectrons <= 0)
      return;
    //auto t6 = std::chrono::high_resolution_clock::now();

    //get electron cluster sizes
    fThisElectronClusterSize = fElectronClusterSize;    
    fNElectronClusters = (size_t)(std::ceil(fNElectrons / fElectronClusterSize));
    if(fNElectronClusters < fMinNumberOfElCluster){
      if(fNElectrons < fMinNumberOfElCluster){
	fNElectronClusters = (size_t)(std::ceil(fNElectrons));
	fThisElectronClusterSize = 1.0;
      }
      else{
	fNElectronClusters = fMinNumberOfElCluster;
	fThisElectronClusterSize = fNElectrons / fMinNumberOfElCluster;
      }
    }
    
    // determine longitudinal & transverse diffusion sigma (cm)
    fLongitudinalDiffSigma = std::sqrt(fTDrift) * fLongitudinalDiffConst;
    fTransverseDiffSigma   = std::sqrt(fTDrift) * fTransverseDiffConst;
    //auto t7 = std::chrono::high_resolution_clock::now();

    //make the diffusion clouds

    //fill electron/energy arrays
    fNElDiff.resize(fNElectronClusters,fThisElectronClusterSize);
    fNEnDiff.resize(fNElectronClusters,fISAlg.EnergyDeposit() * fThisElectronClusterSize/fNElectronClusters);
    fNElDiff.back() = fNElectrons - (fNElectronClusters-1)*fThisElectronClusterSize;
    fNEnDiff.back() = fISAlg.EnergyDeposit() * fNElDiff.back()/fNElectrons;

    //setup position arrays
    fXDiff.resize(fNElectronClusters); fYDiff.resize(fNElectronClusters); fZDiff.resize(fNElectronClusters);
    fTransDiff.resize(2*fNElectronClusters);
    //auto t8 = std::chrono::high_resolution_clock::now();

    // Smear drift times by x position and drift time
    randGauss.fireArray(fNElectronClusters, &fXDiff[0], 0., fLongitudinalDiffSigma);
    
    // Smear the Y,Z position by the transverse diffusion
    randGauss.fireArray(fNElectronClusters, &fYDiff[0], edep.Y()+fPosOffsets[1],fTransverseDiffSigma);
    randGauss.fireArray(fNElectronClusters, &fZDiff[0], edep.Z()+fPosOffsets[2],fTransverseDiffSigma);
    
    //auto t9 = std::chrono::high_resolution_clock::now();

    //now fill the time components
    fTDiff.resize(fNElectronClusters);
    for(size_t i_cl=0; i_cl < fNElectronClusters; ++i_cl)
      fTDiff[i_cl] = fTDrift + fXDiff[i_cl]*fRecipDriftVelocityVolume + edep.T();

    /*
    std::chrono::duration<double, std::milli> tmp;
    tmp = t2 - t1; timer[0] += tmp.count();
    tmp = t3 - t2; timer[1] += tmp.count();
    tmp = t4 - t3; timer[2] += tmp.count();
    tmp = t5 - t4; timer[3] += tmp.count();
    tmp = t6 - t5; timer[4] += tmp.count();
    tmp = t7 - t6; timer[5] += tmp.count();
    tmp = t8 - t7; timer[6] += tmp.count();
    tmp = t9 - t8; timer[7] += tmp.count();
    */
    
  }  

  
  
}// namespace
