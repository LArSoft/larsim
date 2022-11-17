////////////////////////////////////////////////////////////////////////
// Plugin Type: analyzer
// File:        PhotonGen_module.cc
// Description:
// Produce photons at the vertex uniformly distributed in the active volume
// Oct. 20, 2020 by Mu Wei wmu@fnal.gov 2020
////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <cmath>
#include <fstream>
#include <memory>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// ART includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace evgen {
  class PhotonGen : public art::EDProducer {
  public:
    explicit PhotonGen(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);
    void beginRun(art::Run& run);

  private:
    void Sample(simb::MCTruth& truth);

    //Flags
    static const int kUNIF = 0;
    static const int kGAUS = 1;

    //TTree to keep track of where particles have been shot
    bool fFillTree; // Do we want to create a TTree of shot particles?
    TTree* fPhotonGen;
    TLorentzVector fShotPos;
    TLorentzVector fShotMom;
    Int_t fEvID;

    //Scan simulation
    bool fScan; // Do we want to scan x, y, z or just simulate a random point?
    double fPx;
    double fPy; // Fixed coordinate
    double fPz;

    //Parameters used to shoot in distributions
    int fPosDist;   //
    int fTDist;     // Random distributions to use : 1= gauss, 0= uniform
    int fPDist;     //
    double fX;      // central x position of source
    double fY;      // central y position of source
    double fZ;      // central z position of source
    double fT;      // central t position of source
    double fSigmaT; // t width
    double fP;      // central momentm of photon
    double fSigmaP; // mom width;

    // Number of photons per event
    int fN; // number of photons per event

    CLHEP::HepRandomEngine& fEngine;

    //Boundaries of the detector
    double fXmin;
    double fXmax;
    double fYmin;
    double fYmax;
    double fZmin;
    double fZmax;
  };

  //----------------------------------------------------------------
  PhotonGen::PhotonGen(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fFillTree{pset.get<bool>("FillTree")}
    , fScan{pset.get<bool>("Scan")}
    , fPosDist{pset.get<int>("PosDist")}
    , fTDist{pset.get<int>("TDist")}
    , fPDist{pset.get<int>("PDist")}
    , fT{pset.get<double>("T0")}
    , fSigmaT{pset.get<double>("SigmaT")}
    , fP{pset.get<double>("P")}
    , fSigmaP{pset.get<double>("SigmaP")}
    , fN{pset.get<int>("N")}
    , fEngine(art::ServiceHandle<rndm::NuRandomService> {}->createEngine(*this, pset, "Seed"))
  {
    produces<sumdata::RunData, art::InRun>();
    produces<std::vector<simb::MCTruth>>();

    if (fFillTree) {
      art::ServiceHandle<art::TFileService const> tfs;
      fPhotonGen = tfs->make<TTree>("PhGen", "PhGen");
      fPhotonGen->Branch("X", &(fShotPos[0]), "X/D");
      fPhotonGen->Branch("Y", &(fShotPos[1]), "Y/D");
      fPhotonGen->Branch("Z", &(fShotPos[2]), "Z/D");
      fPhotonGen->Branch("T", &(fShotPos[3]), "T/D");
      fPhotonGen->Branch("PX", &(fShotMom[0]), "PX/D");
      fPhotonGen->Branch("PY", &(fShotMom[1]), "PY/D");
      fPhotonGen->Branch("PZ", &(fShotMom[2]), "PZ/D");
      fPhotonGen->Branch("PT", &(fShotMom[3]), "PT/D");
      fPhotonGen->Branch("EventID", &fEvID, "EventID/I");
    }

    if (fScan) {
      fPx = pset.get<double>("Px");
      fPy = pset.get<double>("Py");
      fPz = pset.get<double>("Pz");
      std::cout << "Will generate photons from 3 points." << std::endl;
    }
  }

  //____________________________________________________________________________
  void PhotonGen::beginRun(art::Run& run)
  {
    art::ServiceHandle<geo::Geometry const> geo;
    std::cout << "Number of optical detector: " << int(geo->Cryostat(0).NOpDet()) << std::endl;

    //Boundary of Cryogenic
    double CryoBounds[6];
    geo->CryostatBoundaries(CryoBounds);
    fXmin = CryoBounds[0];
    fXmax = CryoBounds[1];
    fYmin = CryoBounds[2];
    fYmax = CryoBounds[3];
    fZmin = CryoBounds[4];
    fZmax = CryoBounds[5];
    std::cout << "Cryo Boundaries:" << std::endl;
    std::cout << "Xmin: " << fXmin << " Xmax: " << fXmax << " Ymin: " << fYmin << " Ymax: " << fYmax
              << " Zmin: " << fZmin << " Zmax: " << fZmax << std::endl;
    /*
        //Boundary of active volume
        fXmin =  1e9;
        fXmax = -1e9;
        fYmin =  1e9;
        fYmax = -1e9;
        fZmin =  1e9;
        fZmax = -1e9;
        for (size_t i = 0; i < geo->NTPC(); ++ i)
        {
            const geo::TPCGeo &tpc = geo->TPC(i);
            if (fXmin > tpc.MinX()) fXmin = tpc.MinX();
            if (fXmax < tpc.MaxX()) fXmax = tpc.MaxX();
            if (fYmin > tpc.MinY()) fYmin = tpc.MinY();
            if (fYmax < tpc.MaxY()) fYmax = tpc.MaxY();
            if (fZmin > tpc.MinZ()) fZmin = tpc.MinZ();
            if (fZmax < tpc.MaxZ()) fZmax = tpc.MaxZ();
        }
        std::cout << "Active Volume Boundaries:" << std::endl;
        std::cout << "Xmin: " << fXmin << " Xmax: " << fXmax <<" Ymin: " << fYmin << " Ymax: " << fYmax <<" Zmin: " << fZmin << " Zmax: " << fZmax << std::endl;
*/
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));

    return;
  }

  //----------------------------------------------------------------
  void PhotonGen::produce(art::Event& evt)
  {
    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_real_distribution<double> distX(fXmin, fXmax);
    std::uniform_real_distribution<double> distY(fYmin, fYmax);
    std::uniform_real_distribution<double> distZ(fZmin, fZmax);
    std::uniform_real_distribution<double> width(-2.0, 2.0);
    srand((unsigned)time(NULL));

    std::unique_ptr<std::vector<simb::MCTruth>> truthcol(new std::vector<simb::MCTruth>);
    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);

    if (fScan) {
      /*
            int co_index = rand()%3;
            switch(co_index)
            { 
                case 0:
                {
                    fX = distX(rng);
                    fY = fPy + width(rng);
                    fZ = fPz + width(rng);
                    break;
                }
                case 1:
                {
                    fX = fPx + width(rng);
                    fY = distY(rng);
                    fZ = fPz + width(rng);
                    break;
                }
                case 2:
                {
                    fX = fPx + width(rng);
                    fY = fPy + width(rng);
                    fZ = distZ(rng);
                    break;                    
                }
                default:
                {
                    std::cout << "Index Error!" << std::endl;
                    
                    fX = distX(rng);
                    fY = distY(rng);
                    fZ = distZ(rng);
                }
            }
        */
      fX = distX(rng);
      fY = fPy + width(rng);
      fZ = fPz + width(rng);
      Sample(truth);
    }
    else {
      fX = distX(rng);
      fY = distY(rng);
      fZ = distZ(rng);

      Sample(truth);
    }

    truthcol->push_back(truth);
    evt.put(std::move(truthcol));
    return;
  }

  void PhotonGen::Sample(simb::MCTruth& mct)
  {
    std::cout << "Photons Shooting at " << fX << " " << fY << " " << fZ << std::endl;

    CLHEP::RandFlat flat(fEngine);
    CLHEP::RandGaussQ gauss(fEngine);

    for (int j = 0; j != fN; ++j) {
      TVector3 pos;
      pos[0] = fX;
      pos[1] = fY;
      pos[2] = fZ;

      double time;
      if (fTDist == kGAUS) { time = gauss.fire(fT, fSigmaT); }
      else {
        time = fT + fSigmaT * (2.0 * flat.fire() - 1.0);
      }
      fShotPos = TLorentzVector(pos[0], pos[1], pos[2], time);

      //momentum (supplied in eV, convert to GeV)
      double p = fP;
      if (fPDist == kGAUS) { p = gauss.fire(fP, fSigmaP); }
      else {
        p = fP + fSigmaP * (2.0 * flat.fire() - 1.0);
      }
      p /= 1000000000.;

      //angles
      double costh = 2 * flat.fire() - 1;
      double sinth = pow(1 - pow(costh, 2), 0.5);
      double phi = 2 * M_PI * flat.fire();

      //momentum 4-vector
      fShotMom = TLorentzVector(p * sinth * cos(phi), p * sinth * sin(phi), p * costh, p);

      int trackid =
        -1 * (j + 1); // set track id to -i as these are all primary particles and have id <= 0
      int PDG = 0;    //optical photons have PDG 0
      std::string primary("primary");

      simb::MCParticle particle(trackid, PDG, primary);
      particle.AddTrajectoryPoint(fShotPos, fShotMom);

      if (fFillTree) { fPhotonGen->Fill(); }

      mct.Add(particle);
    }

    return;
  }
}

DEFINE_ART_MODULE(evgen::PhotonGen)
