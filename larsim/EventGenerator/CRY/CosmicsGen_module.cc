////////////////////////////////////////////////////////////////////////
/// \file  CosmicsGen_module.cc
/// \brief Generator for cosmic-rays
///
/// Module to produce cosmic ray MC using CRY
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// ROOT includes
#include "TH1.h"
#include "TH2.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nugen/EventGeneratorBase/evgenbase.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nutools/EventGeneratorBase/CRY/CRYHelper.h"

namespace evgen {

  /// A module to check the results from the Monte Carlo generator
  class CosmicsGen : public art::EDProducer {
  public:
    explicit CosmicsGen(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;
    void beginJob() override;
    void beginRun(art::Run& run) override;

    std::vector<double> fbuffbox;

    TH2F* fPhotonAngles;     ///< Photon rate vs angle
    TH2F* fPhotonAnglesLo;   ///< Photon rate vs angle, low momenta
    TH2F* fPhotonAnglesMi;   ///< Photon rate vs angle, middle momenta
    TH2F* fPhotonAnglesHi;   ///< Photon rate vs angle, high momenta
    TH1F* fPhotonCosQ;       ///< Photon rate vs cos(Q)
    TH1F* fPhotonEnergy;     ///< Photon energy (GeV)
    TH1F* fPhotonsPerSample; ///< number of photons in the sampled time window
    TH1F* fPhotonsInCStat;   ///< number of photons in the cryostat during
                             ///< the sampled time window
    TH1F* fPhotonsInTPC;     ///< number of photons in the tpc during
                             ///< the sampled time window

    TH2F* fElectronAngles;     ///< Electron rate vs angle
    TH2F* fElectronAnglesLo;   ///< Electron rate vs angle, low momenta
    TH2F* fElectronAnglesMi;   ///< Electron rate vs angle, middle momenta
    TH2F* fElectronAnglesHi;   ///< Electron rate vs angle, high momenta
    TH1F* fElectronCosQ;       ///< Electron rate vs cos(Q)
    TH1F* fElectronEnergy;     ///< Electron energy (GeV)
    TH1F* fElectronsPerSample; ///< number of electrons in the sampled time window
    TH1F* fElectronsInCStat;   ///< number of electrons in the cryostat during
                               ///< the sampled time window
    TH1F* fElectronsInTPC;     ///< number of electrons in the tpc during
                               ///< the sampled time window

    TH2F* fMuonAngles;               ///< Muon rate vs angle
    TH2F* fMuonAnglesLo;             ///< Muon rate vs angle, low momenta
    TH2F* fMuonAnglesMi;             ///< Muon rate vs angle, middle momenta
    TH2F* fMuonAnglesHi;             ///< Muon rate vs angle, high momenta
    TH1F* fMuonCosQ;                 ///< Muon rate vs cos(Q)
    TH1F* fMuonEnergy;               ///< Muon energy (GeV)
    TH1F* fMuonsPerSample;           ///< number of muons in the sampled time window
    TH1F* fMuonsInCStat;             ///< number of muons in the cryostat during
                                     ///< the sampled time window
    TH1F* fMuonsInTPC;               ///< number of muons in the tpc during
                                     ///< the sampled time window
    CLHEP::HepRandomEngine& fEngine; ///< art-managed random-number engine
    evgb::CRYHelper fCRYHelp;        ///< CRY generator object
  };
}

namespace evgen {

  //____________________________________________________________________________
  CosmicsGen::CosmicsGen(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fbuffbox{pset.get<std::vector<double>>("BufferBox", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0})}
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    , fEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(createEngine(0),
                                                                                 pset,
                                                                                 "Seed"))
    , fCRYHelp{pset, fEngine, art::ServiceHandle<geo::Geometry const>{} -> GetWorldVolumeName()}
  {
    produces<std::vector<simb::MCTruth>>();
    produces<sumdata::RunData, art::InRun>();
  }

  //____________________________________________________________________________
  void CosmicsGen::beginJob()
  {
    art::ServiceHandle<art::TFileService const> tfs;

    fPhotonAngles =
      tfs->make<TH2F>("fPhotonAngles", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fPhotonAnglesLo =
      tfs->make<TH2F>("fPhotonAnglesLo", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fPhotonAnglesMi =
      tfs->make<TH2F>("fPhotonAnglesMi", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fPhotonAnglesHi =
      tfs->make<TH2F>("fPhotonAnglesHi", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fPhotonCosQ = tfs->make<TH1F>("fPhotonCosQ", ";cos#theta;tracks", 50, -1.0, 1.0);
    fPhotonEnergy = tfs->make<TH1F>("fPhotonEnergy", ";E (GeV)", 5000, 0.0, 1000.0);
    fPhotonsPerSample =
      tfs->make<TH1F>("fPhotonsPerSample", ";Number Photons;Samples", 100, 0, 1000);
    fPhotonsInCStat =
      tfs->make<TH1F>("fPhotonsInCryostat", ";Number Photons;Samples", 100, 0, 1000);
    fPhotonsInTPC = tfs->make<TH1F>("fPhotonsInTPC", ";Number Photons;Samples", 100, 0, 1000);

    fElectronAngles =
      tfs->make<TH2F>("fElectronAngles", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fElectronAnglesLo =
      tfs->make<TH2F>("fElectronAnglesLo", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fElectronAnglesMi =
      tfs->make<TH2F>("fElectronAnglesMi", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fElectronAnglesHi =
      tfs->make<TH2F>("fElectronAnglesHi", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fElectronCosQ = tfs->make<TH1F>("fElectronCosQ", ";cos#theta;tracks", 50, -1.0, 1.0);
    fElectronEnergy = tfs->make<TH1F>("fElectronEnergy", ";E (GeV)", 5000, 0.0, 1000.0);
    fElectronsPerSample =
      tfs->make<TH1F>("fElectronsPerSample", ";Number Electrons;Samples", 100, 0, 1000);
    fElectronsInCStat =
      tfs->make<TH1F>("fElectronsInCryotat", ";Number Electrons;Samples", 100, 0, 1000);
    fElectronsInTPC = tfs->make<TH1F>("fElectronsInTPC", ";Number Electrons;Samples", 100, 0, 1000);

    fMuonAngles =
      tfs->make<TH2F>("fMuonAngles", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fMuonAnglesLo =
      tfs->make<TH2F>("fMuonAnglesLo", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fMuonAnglesMi =
      tfs->make<TH2F>("fMuonAnglesMi", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fMuonAnglesHi =
      tfs->make<TH2F>("fMuonAnglesHi", ";#phi;cos#theta", 36, -180.0, 180.0, 50, -1.0, 1.0);
    fMuonCosQ = tfs->make<TH1F>("fMuonCosQ", ";cos#theta;tracks", 50, -1.0, 1.0);
    fMuonEnergy = tfs->make<TH1F>("fMuonEnergy", ";E (GeV)", 5000, 0.0, 1000.0);
    fMuonsPerSample = tfs->make<TH1F>("fMuonsPerSample", ";Number Muons;Samples", 100, 0, 1000);
    fMuonsInCStat = tfs->make<TH1F>("fMuonsInCryostat", ";Number Muons;Samples", 100, 0, 1000);
    fMuonsInTPC = tfs->make<TH1F>("fMuonsInTPC", ";Number Muons;Samples", 100, 0, 1000);
  }

  //____________________________________________________________________________
  void CosmicsGen::beginRun(art::Run& run)
  {
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()), art::fullRun());
  }

  //____________________________________________________________________________
  void CosmicsGen::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<simb::MCTruth>> truthcol(new std::vector<simb::MCTruth>);

    // fill some histograms about this event
    art::ServiceHandle<geo::Geometry const> geom;

    int nCrossCryostat = 0;

    simb::MCTruth truth;

    auto const [surface_y, tpc_length] = std::make_tuple(geom->SurfaceY(), geom->TPC().Length());
    while (nCrossCryostat < 1) {

      simb::MCTruth pretruth;
      truth.SetOrigin(simb::kCosmicRay);
      fCRYHelp.Sample(pretruth, surface_y, tpc_length, 0);

      int numPhotons = 0;
      int numElectrons = 0;
      int numMuons = 0;
      int allPhotons = 0;
      int allElectrons = 0;
      int allMuons = 0;

      // loop over particles in the truth object
      for (int i = 0; i < pretruth.NParticles(); ++i) {
        simb::MCParticle particle = pretruth.GetParticle(i);
        const TLorentzVector& v4 = particle.Position();
        const TLorentzVector& p4 = particle.Momentum();
        double x0[3] = {v4.X(), v4.Y(), v4.Z()};
        double dx[3] = {p4.Px(), p4.Py(), p4.Pz()};

        if (std::abs(particle.PdgCode()) == 13)
          ++allMuons;
        else if (std::abs(particle.PdgCode()) == 22)
          ++allPhotons;
        else if (std::abs(particle.PdgCode()) == 11)
          ++allElectrons;

        TH1F* hCosQ = nullptr;
        TH2F* hAngles = nullptr;
        TH2F* hAnglesLo = nullptr;
        TH2F* hAnglesMi = nullptr;
        TH2F* hAnglesHi = nullptr;
        TH1F* hEnergy = nullptr;
        if (std::abs(particle.PdgCode()) == 13) {
          hCosQ = fMuonCosQ;
          hAngles = fMuonAngles;
          hAnglesLo = fMuonAnglesLo;
          hAnglesMi = fMuonAnglesMi;
          hAnglesHi = fMuonAnglesHi;
          hEnergy = fMuonEnergy;
        }
        else if (std::abs(particle.PdgCode()) == 22) {
          hCosQ = fPhotonCosQ;
          hAngles = fPhotonAngles;
          hAnglesLo = fPhotonAnglesLo;
          hAnglesMi = fPhotonAnglesMi;
          hAnglesHi = fPhotonAnglesHi;
          hEnergy = fPhotonEnergy;
        }
        else if (std::abs(particle.PdgCode()) == 11) {
          hCosQ = fElectronCosQ;
          hAngles = fElectronAngles;
          hAnglesLo = fElectronAnglesLo;
          hAnglesMi = fElectronAnglesMi;
          hAnglesHi = fElectronAnglesHi;
          hEnergy = fElectronEnergy;
        }

        // now check if the particle goes through any cryostat in the detector
        // if so, add it to the truth object.
        for (auto const& cryostat : geom->Iterate<geo::CryostatGeo>()) {

          double bounds[6] = {0.};
          cryostat.Boundaries(bounds);

          //add a buffer box around the cryostat bounds to increase the acceptance
          //(geometrically) at the CRY level to make up for particles we will loose
          //due to multiple scattering effects that pitch in during GEANT4 tracking
          //By default, the buffer box has zero size
          for (unsigned int cb = 0; cb < 6; cb++)
            bounds[cb] = bounds[cb] + fbuffbox[cb];

          //calculate the intersection point with each cryostat surface
          bool intersects_cryo = false;
          for (int bnd = 0; bnd != 6; ++bnd) {
            if (bnd < 2) {
              double p2[3] = {bounds[bnd],
                              x0[1] + (dx[1] / dx[0]) * (bounds[bnd] - x0[0]),
                              x0[2] + (dx[2] / dx[0]) * (bounds[bnd] - x0[0])};
              if (p2[1] >= bounds[2] && p2[1] <= bounds[3] && p2[2] >= bounds[4] &&
                  p2[2] <= bounds[5]) {
                intersects_cryo = true;
                break;
              }
            }
            else if (bnd >= 2 && bnd < 4) {
              double p2[3] = {x0[0] + (dx[0] / dx[1]) * (bounds[bnd] - x0[1]),
                              bounds[bnd],
                              x0[2] + (dx[2] / dx[1]) * (bounds[bnd] - x0[1])};
              if (p2[0] >= bounds[0] && p2[0] <= bounds[1] && p2[2] >= bounds[4] &&
                  p2[2] <= bounds[5]) {
                intersects_cryo = true;
                break;
              }
            }
            else if (bnd >= 4) {
              double p2[3] = {x0[0] + (dx[0] / dx[2]) * (bounds[bnd] - x0[2]),
                              x0[1] + (dx[1] / dx[2]) * (bounds[bnd] - x0[2]),
                              bounds[bnd]};
              if (p2[0] >= bounds[0] && p2[0] <= bounds[1] && p2[1] >= bounds[2] &&
                  p2[1] <= bounds[3]) {
                intersects_cryo = true;
                break;
              }
            }
          }

          if (intersects_cryo) {
            truth.Add(particle);

            if (std::abs(particle.PdgCode()) == 13)
              ++numMuons;
            else if (std::abs(particle.PdgCode()) == 22)
              ++numPhotons;
            else if (std::abs(particle.PdgCode()) == 11)
              ++numElectrons;

            if (hCosQ != 0) {
              double cosq = -p4.Py() / p4.P();
              double phi = std::atan2(p4.Pz(), p4.Px());
              phi *= 180 / M_PI;
              hCosQ->Fill(cosq);
              hAngles->Fill(phi, cosq);
              if (p4.E() < 1.0)
                hAnglesLo->Fill(phi, cosq);
              else if (p4.E() < 10.0)
                hAnglesMi->Fill(phi, cosq);
              else
                hAnglesHi->Fill(phi, cosq);
              hEnergy->Fill(p4.E());
            }      //end if there is a cos(theta) histogram
            break; //leave loop over cryostats to avoid adding particle multiple times
          }        // end if particle goes into a cryostat
        }          // end loop over cryostats in the detector

      } // loop on particles

      nCrossCryostat = truth.NParticles();

      fPhotonsPerSample->Fill(allPhotons);
      fElectronsPerSample->Fill(allElectrons);
      fMuonsPerSample->Fill(allMuons);

      fPhotonsInCStat->Fill(numPhotons);
      fElectronsInCStat->Fill(numElectrons);
      fMuonsInCStat->Fill(numMuons);
    }

    truthcol->push_back(truth);
    evt.put(std::move(truthcol));
  } // end produce

} // end namespace

DEFINE_ART_MODULE(evgen::CosmicsGen)
