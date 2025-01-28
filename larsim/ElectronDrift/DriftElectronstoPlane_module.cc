/**
 * @file SimDriftElectrons_module.cxx
 *
 * @brief Transports energy depositions in the LAr TPC to the TPC
 * channels.
 *
 * author:
 * This module was prepared by William Seligman (me), based on code
 * that had been in
 * `LArG4::LArVoxelReadout::DriftIonizationElectrons`. However, though
 * I wrote the original LArVoxelReadout code, I have no idea who added
 * DriftIonizationElectrons. I probably will not be able to answer any
 * questions about how this code works.
 *
 * This module acts on sim::SimEnergyDeposit, the single energy
 * depositions from the detector simulation (LArG4), and simulates the
 * transport of the ensuing ionization electrons to the readout
 * channels:
 *
 * 1. the number of ionisation electrons is read from the current
 *   `larg4::IonizationAndScintillation` instance
 * 2. space charge displacement is optionally applied
 * 3. lifetime correction is applied
 * 4. charge is split in small electron clusters
 * 5. each cluster is subject to longitudinal and transverse diffusion
 * 6. each cluster is assigned to one TPC channel for each wire plane
 * 7. optionally, charge is forced to stay on the planes; otherwise charge
 *    drifting outside the plane is lost
 *
 * For each energy deposition, entries on the appropriate
 * `sim::SimChannel` are added, with the information of the position
 * where the energy deposit happened (in global coordinates,
 * centimeters), the ID of the track in the detector simulation which
 * produced the deposition, and the quantized time of arrival to the
 * channel (in global TDC tick units).  At most one entry is added for
 * each electron cluster, but entries from the same energy deposit can
 * be compacted if falling on the same TDC tick.
 *
 * Options
 * --------
 *
 * A few optional behaviours are supported:
 *
 * * lead off-plane charge to the planes: regulated by
 *   `RecoverOffPlaneDeposit()`, if charge which reaches a wire plane
 *   is actually off it by less than the chosen margin, it's accounted for by
 *   that plane; by default the margin is 0 and all the charge off the plane
 *   is lost (with a warning)
 *
 * Update:
 * Christoph Alt, September 2018 (christoph.alt@cern.ch)
 * Break hardcoded charge drift in x to support charge drift in y and z.
 */

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/Simulation/SimDriftedElectronCluster.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/ElectronDrift/ISCalculationSeparate.h"
#include "larsim/Utils/SCEOffsetBounds.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"

// External libraries
#include "CLHEP/Random/RandGauss.h"
#include "TMath.h"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace detsim {

  // Base class for creation of raw signals on wires.
  class DriftElectronstoPlane : public art::EDProducer {
  public:
    explicit DriftElectronstoPlane(fhicl::ParameterSet const& pset);

    // Methods that that are available for a module derived from
    // art::EDProducer.
    void produce(art::Event& evt) override;
    void beginJob() override;

  private:
    // The label of the module that created the sim::SimEnergyDeposit
    // objects (as of Oct-2017, this is probably "largeant").
    art::InputTag fSimModuleLabel;

    CLHEP::RandGauss fRandGauss;

    bool fStoreDriftedElectronClusters;
    double fLongitudinalDiffusion;
    double fTransverseDiffusion;
    double fElectronClusterSize;
    int fMinNumberOfElCluster;
    double fGeVToElectrons;
    double fRecombA;
    double fRecombk;
    double fModBoxA;
    double fModBoxB;
    bool fUseModBoxRecomb;

    double fElectronLifetime;
    double fLifetimeCorr_const;
    double fLDiff_const;
    double fTDiff_const;
    double fRecipDriftVel[3];

    // Save the number of cryostats, and the number of TPCs within
    // each cryostat.
    size_t fNCryostats;
    std::vector<size_t> fNTPCs;

    // Per-cluster information.
    std::vector<double> fLongDiff;
    std::vector<double> fTransDiff1;
    std::vector<double> fTransDiff2;
    std::vector<double> fnElDiff;
    std::vector<double> fnEnDiff;

    double fDriftClusterPos[3];

    art::ServiceHandle<geo::Geometry const> fGeometry; ///< Handle to the Geometry service

    //IS calculationg
    ISCalculationSeparate fISAlg;

  }; // class DriftElectronstoPlane

  //-------------------------------------------------
  DriftElectronstoPlane::DriftElectronstoPlane(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fSimModuleLabel{pset.get<art::InputTag>("SimulationLabel")}
    // create a default random engine; obtain the random seed from
    // NuRandomService, unless overridden in configuration with key
    // "Seed"
    , fRandGauss{art::ServiceHandle<rndm::NuRandomService>{}
                 -> registerAndSeedEngine(createEngine(0), pset, "Seed")}
    , fStoreDriftedElectronClusters{pset.get<bool>("StoreDriftedElectronClusters", true)}
    , fLongitudinalDiffusion{pset.get<double>("LongitudinalDiffusion", 6.2e-9)}
    , fTransverseDiffusion{pset.get<double>("TransverseDiffusion", 16.3e-9)}
    , fElectronClusterSize{pset.get<double>("ElectronClusterSize", 600.0)}
    , fMinNumberOfElCluster{pset.get<int>("MinNumberOfElCluster", 0)}
    , fGeVToElectrons{pset.get<double>("GeVToElectrons", 4.237e+07)}
    , fRecombA{pset.get<double>("RecombA", 0.800)}
    , fRecombk{pset.get<double>("Recombk", 0.0486)}
    , fModBoxA{pset.get<double>("ModBoxA", 0.930)}
    , fModBoxB{pset.get<double>("ModBoxB", 0.212)}
    , fUseModBoxRecomb{pset.get<bool>("UseModBoxRecomb", true)}
    , fISAlg{pset}
  {
    if (fStoreDriftedElectronClusters) { produces<std::vector<sim::SimDriftedElectronCluster>>(); }
  }

  //-------------------------------------------------
  void DriftElectronstoPlane::beginJob()
  {
    // Define the physical constants we'll use.

    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    fElectronLifetime =
      detProp
        .ElectronLifetime(); // Electron lifetime as returned by the DetectorProperties service assumed to be in us;
    for (int i = 0; i < 3; ++i) {
      double driftVelocity =
        detProp.DriftVelocity(detProp.Efield(i),
                              detProp.Temperature()) *
        1.e-3; //  Drift velocity as returned by the DetectorProperties service assumed to be in cm/us. Multiply by 1.e-3 to convert into LArSoft standard velocity units, cm/ns;

      fRecipDriftVel[i] = 1. / driftVelocity;
    }
    MF_LOG_DEBUG("DriftElectronstoPlane")
      << " e lifetime (ns): " << fElectronLifetime
      << "\n Temperature (K): " << detProp.Temperature()
      << "\n Drift velocity (cm/ns): " << 1. / fRecipDriftVel[0] << " " << 1. / fRecipDriftVel[1]
      << " " << 1. / fRecipDriftVel[2];

    // Opposite of lifetime. Convert from us to standard LArSoft time units, ns;
    fLifetimeCorr_const = -1000. * fElectronLifetime;
    fLDiff_const = std::sqrt(2. * fLongitudinalDiffusion);
    fTDiff_const = std::sqrt(2. * fTransverseDiffusion);

    // For this detector's geometry, save the number of cryostats and the number of TPCs
    // within each cryostat.
    fNCryostats = fGeometry->Ncryostats();
    fNTPCs.resize(fNCryostats);
    for (size_t n = 0; n < fNCryostats; ++n)
      fNTPCs[n] = fGeometry->NTPC(geo::CryostatID(n));
  }

  //-------------------------------------------------
  void DriftElectronstoPlane::produce(art::Event& event)
  {
    // Fetch the SimEnergyDeposit objects for this event.
    typedef art::Handle<std::vector<sim::SimEnergyDeposit>> energyDepositHandle_t;
    energyDepositHandle_t energyDepositHandle;
    // If there aren't any energy deposits for this event, don't panic. It's possible
    // someone is doing a study with events outside the TPC, or where there are only
    // non-ionizing particles, or something like that.
    if (!event.getByLabel(fSimModuleLabel, energyDepositHandle)) return;
    // Container for the SimDriftedElectronCluster objects
    std::unique_ptr<std::vector<sim::SimDriftedElectronCluster>>
      SimDriftedElectronClusterCollection(new std::vector<sim::SimDriftedElectronCluster>);

    // We're going through the input vector by index, rather than by iterator, because we
    // need the index number to compute the associations near the end of this method.
    auto const& energyDeposits = *energyDepositHandle;
    auto energyDepositsSize = energyDeposits.size();

    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);

    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();
    // For each energy deposit in this event
    for (size_t edIndex = 0; edIndex < energyDepositsSize; ++edIndex) {
      auto const& energyDeposit = energyDeposits[edIndex];

      // "xyz" is the position of the energy deposit in world coordinates. Note that the
      // units of distance in sim::SimEnergyDeposit are supposed to be cm.
      auto const mp = energyDeposit.MidPoint();
      double const xyz[3] = {mp.X(), mp.Y(), mp.Z()};

      // From the position in world coordinates, determine the cryostat and tpc. If
      // somehow the step is outside a tpc (e.g., cosmic rays in rock) just move on to the
      // next one.
      geo::TPCGeo const& tpcGeo = fGeometry->TPC();
      geo::TPCID const& tpcid = tpcGeo.ID();
      unsigned const nplanes = wireReadoutGeom.Nplanes(tpcid);

      // The drift direction can be either in the positive or negative direction in any
      // coordinate x, y or z.

      // Define charge drift direction: driftcoordinate (x, y or z) and driftsign
      // (positive or negative). Also define coordinates perpendicular to drift direction.
      auto const [axis, sign] = tpcGeo.DriftAxisWithSign();

      // For axis == geo::Coordinate::X
      int driftcoordinate = 0;
      int transversecoordinate1 = 1;
      int transversecoordinate2 = 2;
      if (axis == geo::Coordinate::Y) {
        transversecoordinate1 = 0;
        driftcoordinate = 1;
        transversecoordinate2 = 2;
      }
      else if (axis == geo::Coordinate::Z) {
        transversecoordinate1 = 0;
        transversecoordinate2 = 1;
        driftcoordinate = 2;
      }

      int const driftsign = to_int(sign); //1: +x, +y or +z, -1: -x, -y or -z

      //Check for charge deposits behind charge readout planes
      auto const plane_center = wireReadoutGeom.FirstPlane(tpcid).GetCenter();
      auto const plane_center_coord = geo::vect::coord(plane_center, driftcoordinate);
      if (driftsign == 1 && plane_center_coord < xyz[driftcoordinate]) continue;
      if (driftsign == -1 && plane_center_coord > xyz[driftcoordinate]) continue;

      /// \todo think about effects of drift between planes.
      // Center of plane is also returned in cm units
      double DriftDistance = std::abs(xyz[driftcoordinate] - plane_center_coord);

      // Space-charge effect (SCE): Get SCE {x,y,z} offsets for particular location in TPC
      geo::Vector_t posOffsets{0.0, 0.0, 0.0};
      double posOffsetxyz[3] = {
        0.0, 0.0, 0.0}; //need this array for the driftcoordinate and transversecoordinates
      auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
      if (SCE->EnableSimSpatialSCE() == true) {
        posOffsets = SCE->GetPosOffsets(mp);
        if (larsim::Utils::SCE::out_of_bounds(posOffsets)) continue;
        posOffsetxyz[0] = posOffsets.X();
        posOffsetxyz[1] = posOffsets.Y();
        posOffsetxyz[2] = posOffsets.Z();
      }

      double avegagetransversePos1 = 0.;
      double avegagetransversePos2 = 0.;

      DriftDistance += -1. * posOffsetxyz[driftcoordinate];
      avegagetransversePos1 = xyz[transversecoordinate1] + posOffsetxyz[transversecoordinate1];
      avegagetransversePos2 = xyz[transversecoordinate2] + posOffsetxyz[transversecoordinate2];

      // Space charge distortion could push the energy deposit beyond the wire plane (see
      // issue #15131). Given that we don't have any subtlety in the simulation of this
      // region, bringing the deposit exactly on the plane should be enough for the time
      // being.
      if (DriftDistance < 0.) DriftDistance = 0.;

      // Drift time in ns
      double TDrift = DriftDistance * fRecipDriftVel[0];

      if (nplanes == 2 &&
          driftcoordinate ==
            0) { // special case for ArgoNeuT (Nplanes = 2 and drift direction = x): plane 0
                 // is the second wire plane
        auto const pitch = wireReadoutGeom.PlanePitch(tpcGeo.ID(), 0, 1);
        TDrift = ((DriftDistance - pitch) * fRecipDriftVel[0] + pitch * fRecipDriftVel[1]);
      }
      const int nIonizedElectrons =
        fISAlg.CalculateIonizationAndScintillation(detProp, energyDeposit).numElectrons;

      const double lifetimecorrection = TMath::Exp(TDrift / fLifetimeCorr_const);
      const double energy = energyDeposit.Energy();

      // if we have no electrons (too small energy or too large recombination) we are done
      // already here
      if (nIonizedElectrons <= 0) {
        MF_LOG_DEBUG("DriftElectronstoPlane")
          << "step " // << energyDeposit << "\n"
          << "No electrons drifted to readout, " << energy << " MeV lost.";
        continue;
      }

      // includes the effect of lifetime: lifetimecorrection = exp[-tdrift/tau]
      const double nElectrons = nIonizedElectrons * lifetimecorrection;

      // Longitudinal & transverse diffusion sigma (cm)
      double SqrtT = std::sqrt(TDrift);
      double LDiffSig = SqrtT * fLDiff_const;
      double TDiffSig = SqrtT * fTDiff_const;
      double electronclsize = fElectronClusterSize;

      // Number of electron clusters.
      int nClus = (int)std::ceil(nElectrons / electronclsize);
      if (nClus < fMinNumberOfElCluster) {
        electronclsize = nElectrons / fMinNumberOfElCluster;
        if (electronclsize < 1.0) { electronclsize = 1.0; }
        nClus = (int)std::ceil(nElectrons / electronclsize);
      }

      // Empty and resize the electron-cluster vectors.
      fLongDiff.clear();
      fTransDiff1.clear();
      fTransDiff2.clear();
      fnElDiff.clear();
      fnEnDiff.clear();
      fLongDiff.resize(nClus);
      fTransDiff1.resize(nClus);
      fTransDiff2.resize(nClus);
      fnElDiff.resize(nClus, electronclsize);
      fnEnDiff.resize(nClus);

      // fix the number of electrons in the last cluster, that has a smaller size
      fnElDiff.back() = nElectrons - (nClus - 1) * electronclsize;

      for (size_t xx = 0; xx < fnElDiff.size(); ++xx) {
        if (nElectrons > 0)
          fnEnDiff[xx] = energy / nElectrons * fnElDiff[xx];
        else
          fnEnDiff[xx] = 0.;
      }

      // Smear drift times by longitudinal diffusion
      if (LDiffSig > 0.0)
        fRandGauss.fireArray(nClus, &fLongDiff[0], 0., LDiffSig);
      else
        fLongDiff.assign(nClus, 0.0);

      if (TDiffSig > 0.0) {
        // Smear the coordinates in plane perpendicular to drift direction by the transverse diffusion
        fRandGauss.fireArray(nClus, &fTransDiff1[0], avegagetransversePos1, TDiffSig);
        fRandGauss.fireArray(nClus, &fTransDiff2[0], avegagetransversePos2, TDiffSig);
      }
      else {
        fTransDiff1.assign(nClus, avegagetransversePos1);
        fTransDiff2.assign(nClus, avegagetransversePos2);
      }

      // make a collection of electrons for each plane
      for (auto const& plane : wireReadoutGeom.Iterate<geo::PlaneGeo>(tpcGeo.ID())) {
        auto const plane_center = plane.GetCenter();
        fDriftClusterPos[driftcoordinate] = geo::vect::coord(plane_center, driftcoordinate);

        // Drift nClus electron clusters to the induction plane
        for (int k = 0; k < nClus; ++k) {

          // Correct drift time for longitudinal diffusion and plane
          double TDiff = TDrift + fLongDiff[k] * fRecipDriftVel[0];

          // Take into account different Efields between planes.  Also take into account
          // special case for ArgoNeuT (Nplanes = 2 and drift direction = x): plane 0 is
          // the second wire plane
          for (unsigned int ip = 0; ip < plane.ID().Plane; ++ip) {
            auto const plane_center = wireReadoutGeom.Plane({tpcGeo.ID(), ip}).GetCenter();
            auto const next_plane_center = wireReadoutGeom.Plane({tpcGeo.ID(), ip + 1}).GetCenter();
            TDiff += (geo::vect::coord(next_plane_center, driftcoordinate) -
                      geo::vect::coord(plane_center, driftcoordinate)) *
                     fRecipDriftVel[(nplanes == 2 && driftcoordinate == 0) ? ip + 2 : ip + 1];
          }

          fDriftClusterPos[transversecoordinate1] = fTransDiff1[k];
          fDriftClusterPos[transversecoordinate2] = fTransDiff2[k];
          auto const simTime = energyDeposit.Time();
          /// \todo think about effects of drift between planes
          SimDriftedElectronClusterCollection->emplace_back(
            fnElDiff[k],
            TDiff + simTime,                      // timing
            geo::Point_t{mp.X(), mp.Y(), mp.Z()}, // mean position of the deposited energy
            geo::Point_t{fDriftClusterPos[0],
                         fDriftClusterPos[1],
                         fDriftClusterPos[2]}, // final position of the drifted cluster
            geo::Point_t{
              LDiffSig, TDiffSig, TDiffSig}, // Longitudinal (X) and transverse (Y,Z) diffusion
            fnEnDiff[k],                     //deposited energy that originated this cluster
            energyDeposit.TrackID());

        } // end loop over clusters
      }   // end loop over planes
    }     // for each sim::SimEnergyDeposit

    if (fStoreDriftedElectronClusters) event.put(std::move(SimDriftedElectronClusterCollection));
  }
} // namespace detsim

DEFINE_ART_MODULE(detsim::DriftElectronstoPlane)
