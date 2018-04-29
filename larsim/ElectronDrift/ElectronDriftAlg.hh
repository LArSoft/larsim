////////////////////////////////////////////////////////////////////////
/// \file  ElectronDrift.hh
/// \brief Algorithm class for modeling electron drift, given an
///        energy deposition in space.
///
/// \version $Id:  $
/// \author  wketchum@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_ELECTRONDRIFT_H
#define LARG4_ELECTRONDRIFT_H

#include <unordered_map>

// forward declaration
namespace detinfo { class LArProperties; class DetectorProperties; }
namespace sim { class SimEnergyDeposit; class LArG4Parameters; }
namespace spacecharge { class SpaceCharge; }
namespace geo { class GeometryCore; }

#include "larsim/IonizationScintillation/ISCalculationSeparate.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "CLHEP/Random/RandomEngine.h"


namespace larg4 {
  
  class ElectronDriftAlg{
    
  public:
    
   ElectronDriftAlg();
    virtual ~ElectronDriftAlg();
    
    void   Initialize(const detinfo::LArProperties* larp,   
		      const detinfo::DetectorProperties* detp,
		      const sim::LArG4Parameters* lgp,
		      const spacecharge::SpaceCharge* sce,
		      const geo::GeometryCore* geo,
		      CLHEP::HepRandomEngine* engine);

    size_t const& NElectronClusters() const { return fNElectronClusters; }
    std::vector<double> const& ElectronClustersY()  const { return fYDiff; }
    std::vector<double> const& ElectronClustersZ()  const { return fZDiff; }
    std::vector<double> const& ElectronClustersX()  const { return fXDiff; }
    std::vector<double> const& ElectronClustersT()  const { return fTDiff; }
    std::vector<double> const& ElectronClustersEl() const { return fNElDiff; }
    std::vector<double> const& ElectronClustersEn() const { return fNEnDiff; }

    std::vector<double> const& GetTimer() const { return timer; }
    
    void PrecalculateSC(std::vector<sim::SimEnergyDeposit> const&);
    void DriftElectrons(sim::SimEnergyDeposit const& edep);

    
  private:
    
    const spacecharge::SpaceCharge*  fSCE;
    const geo::GeometryCore*         fGeo;
    
    double fElectronLifetime;
    double fDriftVelocityVolume;
    double fRecipDriftVelocityVolume;
    double fLongitudinalDiffConst;
    double fLongitudinalDiffSigma;
    double fTransverseDiffConst;
    double fTransverseDiffSigma;
    double fElectronClusterSize;
    double fMinNumberOfElCluster;

    
    CLHEP::HepRandomEngine* fEngine;
    
    double fXYZ[3];
    geo::TPCID fTPC_id;
    double fXDrift;
    double fTDrift;
    std::vector<double> fPosOffsets;
    
    ISCalculationSeparate fISAlg;
    double fNElectrons;
    size_t fNElectronClusters;
    double fThisElectronClusterSize;

    std::vector<double> fXDiff;
    std::vector<double> fYDiff;
    std::vector<double> fZDiff;
    std::vector<double> fTransDiff;
    std::vector<double> fTDiff;
    std::vector<double> fNElDiff;
    std::vector<double> fNEnDiff;

    std::vector<double> timer;

    unsigned long N_VOXELS_X;
    unsigned long N_VOXELS_Y;
    unsigned long N_VOXELS_Z;
    std::unordered_map< unsigned long, std::vector< std::vector<double> > > fSCCalcMap;

    unsigned long GetSCMapIndex(float,float,float);
    
  };
}
#endif // LARG4_ELECTRONDRIFT_H

