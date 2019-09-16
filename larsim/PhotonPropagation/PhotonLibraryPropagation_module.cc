/**
 * @file larsim/PhotonPropagation/PhotonLibraryPropagation_module.cc
 * @brief Provides `phot:PhotonLibraryPropagation` module.
 */

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include <cmath>
#include <memory>

using namespace std;

namespace {

double single_exp(double t, double tau2)
{
         return exp((-1.0 * t) / tau2) / tau2;
}

double bi_exp(double t, double tau1, double tau2)
{
         return (((exp((-1.0 * t) / tau2) * (1.0 - exp((-1.0 * t) / tau1))) / tau2) / tau2) * (tau1 + tau2);
}


// Returns the time within the time distribution of the scintillation process, when the photon was created.
// Scintillation light has an exponential decay which here is given by the decay time, tau2,
// and an exponential increase, which here is given by the rise time, tau1.
// randflatscinttime is passed to use the saved seed from the RandomNumberSaver in order to be able to reproduce the same results.
double GetScintTime(double tau1, double tau2, CLHEP::RandFlat& randflatscinttime)
{
  // tau1: rise time (originally defaulted to -1) and tau2: decay time
  //ran1, ran2 = random numbers for the algorithm
  if ((tau1 == 0.0) || (tau1 == -1.0)) {
    return -tau2 * log(randflatscinttime());
  }
  while (1) {
    auto ran1 = randflatscinttime();
    auto ran2 = randflatscinttime();
    auto d = (tau1 + tau2) / tau2;
    auto t = -tau2 * log(1 - ran1);
    auto g = d * single_exp(t, tau2);
    if (ran2 <= bi_exp(t, tau1, tau2) / g) {
      return t;
    }
  }
}

double GetScintYield(sim::SimEnergyDeposit const& edep, detinfo::LArProperties const& larp)
{
  if (!larp.ScintByParticleType()) {
    return larp.ScintYieldRatio();
  }
  switch (edep.PdgCode()) {
    case 2212:
      return larp.ProtonScintYieldRatio();
    case 13:
    case -13:
      return larp.MuonScintYieldRatio();
    case 211:
    case -211:
      return larp.PionScintYieldRatio();
    case 321:
    case -321:
      return larp.KaonScintYieldRatio();
    case 1000020040:
      return larp.AlphaScintYieldRatio();
    case 11:
    case -11:
    case 22:
      return larp.ElectronScintYieldRatio();
    default:
      return larp.ElectronScintYieldRatio();
  }
}

} // unnamed namespace

namespace phot {

/**
 * @brief Fast simulation of propagating the photons created from SimEnergyDeposits.
 *
 * This module does a fast simulation of propagating the photons created from SimEnergyDeposits,
 * which is the Geant4 output after each step, to each of the optical detectors.
 * This simulation is done using the PhotonLibrary,
 * which stores the visibilities of each optical channel with respect to each optical voxel in the TPC volume,
 * to avoid propagating single photons using Geant4.
 * At the end of this module a collection of the propagated photons either as `sim::SimPhotonsLite` or as `sim::SimPhotons` is placed into the art event.
 *
 * Keep in mind that at this stage the LArG4 main module is not capable of running the full optical simulation,
 * because the necessary code has not yet been written.
 *
 * In the future when the PhotonLibrary has the propagation time included,
 * it could be possible to enhance `sim::SimPhotons` and `sim::SimPhotonsLite` to contain the propagation time.
 * At this point the time recorded for the photon is the creation time of the photon.
 *
 * The steps this module takes are:
 *   - to take `sim::SimEnergyDeposits` produced by larg4Main,
 *   - use Ionisation and Scintillation to calculate the amount of scintillated photons,
 *   - use the PhotonLibrary (visibilities) to determine the amount of visible photons at each optical channel,
 *     - visible photons = the amount of scintillated photons calculated times the visibility
 *       at the middle of the Geant4 step for a given optical channel.
 *   - and if `sim::SimPhotonsLite` produced:
 *         - since a `sim::SimPhotonsLite` only keeps a set of times with the number of photons produced at each time
 *           for a given OpChannel number:
 *              - for each time (as an integer in [ns]) photons are produced along the Geant4 step
 *                (taking into account the rise time and decay time of the fast and slow components of the scintillation process),
 *              - count the amount of photons emitted at that time.
 *          - the total amount of visible photons produced during the current Geant4 step equals the sum of counts for each time.
 *          - the total amount of visible photons produced during the current Geant4 step
 *            is determined by throwing a random number from a Poisson distribution
 *            with a mean of the amount of visible photons calculated above.
 *
 *   - and if `sim::SimPhotons` produced:
 *         - since a `sim::SimPhotons` keeps a collection of photons (`sim::OnePhoton`) for a given OpChannel number:
 *           - each single photon produced by this algorithm is just a copy containing the same information about:
 *             - energy (set to a constant value = 9.7e-6, which is equivalent to a wavelength of 128 nm,
 *               it should actually be 126.6 nm!!),
 *             - initial position,
 *             - time (as an integer in [ns]) the photon is produced along the Geant4 Step
 *               (taking into account the rise time and decay time of the fast and slow components of the scintillation process),
 *           - the total amount of photon copies produced during the current Geant4 step
 *             is determined by throwing a random number from a Poisson distribution
 *             with a mean of the amount of visible photons calculated above.
 *
 * This module should only be run for the fast optical simulation even though it can create `sim::SimPhotonsLite` and `sim::SimPhotons` as data products.
 * If there is need to create `sim::SimPhotons`, there are some considerations you must be aware of.
 * Since the amount of `sim::SimPhotons` produced even at low energies and in small geometries quickly exceeds the memory capacity of the job,
 * right now it is actually impossible to produce `sim::SimPhotons` for any realistic geometry.
 * A possible way around the problem is to implement a scaling of the produced `sim::SimPhotons`, to only produce a fraction of them.
 */
class PhotonLibraryPropagation : public art::EDProducer {

private:

  double fRiseTimeFast;
  double fRiseTimeSlow;
  bool fDoSlowComponent;
  vector<art::InputTag> fEDepTags;
  larg4::ISCalcSeparate fISAlg;
  CLHEP::HepRandomEngine& fPhotonEngine;
  CLHEP::HepRandomEngine& fScintTimeEngine;

  void produce(art::Event&) override;

public:

  explicit PhotonLibraryPropagation(fhicl::ParameterSet const&);
  PhotonLibraryPropagation(PhotonLibraryPropagation const&) = delete;
  PhotonLibraryPropagation(PhotonLibraryPropagation&&) = delete;
  PhotonLibraryPropagation& operator=(PhotonLibraryPropagation const&) = delete;
  PhotonLibraryPropagation& operator=(PhotonLibraryPropagation&&) = delete;
};

PhotonLibraryPropagation::PhotonLibraryPropagation(fhicl::ParameterSet const& p)
  : art::EDProducer{p}
  , fRiseTimeFast{p.get<double>("RiseTimeFast", 0.0)}
  , fRiseTimeSlow{p.get<double>("RiseTimeSlow", 0.0)}
  , fDoSlowComponent{p.get<bool>("DoSlowComponent")}
  , fEDepTags{p.get<vector<art::InputTag>>("EDepModuleLabels")}
  , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "photon",    p, "SeedPhoton"))
  , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "scinttime", p, "SeedScintTime"))
{
  if (art::ServiceHandle<sim::LArG4Parameters const>{}->UseLitePhotons()) {
    produces<vector<sim::SimPhotonsLite>>();
  }
  else {
    produces<vector<sim::SimPhotons>>();
  }
}

void PhotonLibraryPropagation::produce(art::Event& e)
{
  art::ServiceHandle<PhotonVisibilityService const> pvs;
  art::ServiceHandle<sim::LArG4Parameters const> lgp;
  auto const* larp = lar::providerFrom<detinfo::LArPropertiesService>();
  CLHEP::RandPoissonQ randpoisphot{fPhotonEngine};
  CLHEP::RandFlat randflatscinttime{fScintTimeEngine};
  auto const nOpChannels = pvs->NOpChannels();
  fISAlg.Initialize(larp, lar::providerFrom<detinfo::DetectorPropertiesService>(), &*lgp, lar::providerFrom<spacecharge::SpaceChargeService>());
  unique_ptr<vector<sim::SimPhotons>> photCol{new vector<sim::SimPhotons>{}};
  auto& photonCollection{*photCol};
  photonCollection.resize(nOpChannels);
  unique_ptr<vector<sim::SimPhotonsLite>> photLiteCol{new vector<sim::SimPhotonsLite>{}};
  auto& photonLiteCollection{*photLiteCol};
  photonLiteCollection.resize(nOpChannels);
  for (unsigned int i = 0; i < nOpChannels; ++i) {
    photonLiteCollection[i].OpChannel = i;
    photonCollection[i].fOpChannel = i;
  }
  vector<vector<sim::SimEnergyDeposit> const*> edep_vecs;
  for (auto label : fEDepTags) {
    auto const& edep_handle = e.getValidHandle<vector<sim::SimEnergyDeposit>>(label);
    edep_vecs.push_back(edep_handle);
  }
  for (auto const& edeps : edep_vecs) { //loop over modules
    for (auto const& edep : *edeps) { //loop over energy deposits: one per step
      //int count_onePhot =0; // unused
      auto const& p = edep.MidPoint();
      auto const& Visibilities = pvs->GetAllVisibilities(p);
      if (!Visibilities) {
        throw cet::exception("PhotonLibraryPropagation")
          << "There is no entry in the PhotonLibrary for this position in space. "
             "Position: " << edep.MidPoint();
      }
      fISAlg.Reset();
      fISAlg.CalculateIonizationAndScintillation(edep);
      //total amount of scintillation photons
      double nphot = static_cast<int>(fISAlg.NumberScintillationPhotons());
      //amount of scintillated photons created via the fast scintillation process
      double nphot_fast = static_cast<int>(GetScintYield(edep, *larp) * nphot);
      //amount of scintillated photons created via the slow scintillation process
      double nphot_slow = nphot - nphot_fast;
      for (unsigned int channel = 0; channel < nOpChannels; ++channel) {
        auto visibleFraction = Visibilities[channel];
        if (visibleFraction == 0.0) {
          // Voxel is not visible at this optical channel, skip doing anything for this channel.
          continue;
        }
        if (lgp->UseLitePhotons()) {
          if (nphot_fast > 0) {
            //throwing a random number from a poisson distribution with a mean of the amount of photons visible at this channel
            auto n = static_cast<int>(randpoisphot.fire(nphot_fast * visibleFraction));
            for (long i = 0; i < n; ++i) {
              //calculates the time at which the photon was produced
              auto time = static_cast<int>(edep.T0() + GetScintTime(fRiseTimeFast, larp->ScintFastTimeConst(), randflatscinttime));
              ++photonLiteCollection[channel].DetectedPhotons[time];
            }
          }
          if ((nphot_slow > 0) && fDoSlowComponent) {
            //throwing a random number from a poisson distribution with a mean of the amount of photons visible at this channel
            auto n = randpoisphot.fire(nphot_slow * visibleFraction);
            for (long i = 0; i < n; ++i) {
              //calculates the time at which the photon was produced
              auto time = static_cast<int>(edep.T0() + GetScintTime(fRiseTimeSlow, larp->ScintSlowTimeConst(), randflatscinttime));
              ++photonLiteCollection[channel].DetectedPhotons[time];
            }
          }
        }
        else {
          sim::OnePhoton photon;
          photon.SetInSD = false;
          photon.InitialPosition = TVector3{edep.X(), edep.Y(), edep.Z()};
          photon.Energy = 9.7e-6;
          if (nphot_fast > 0) {
            //throwing a random number from a poisson distribution with a mean of the amount of photons visible at this channel
            auto n = randpoisphot.fire(nphot_fast * visibleFraction);
            if (n > 0) {
              //calculates the time at which the photon was produced
              photon.Time = edep.T0() + GetScintTime(fRiseTimeFast, larp->ScintFastTimeConst(), randflatscinttime);
              // add n copies of sim::OnePhoton photon to the photon collection for a given OpChannel
              photonCollection[channel].insert(photonCollection[channel].end(), n, photon);
            }
          }
          if ((nphot_slow > 0) && fDoSlowComponent) {
            //throwing a random number from a poisson distribution with a mean of the amount of photons visible at this channel
            auto n = randpoisphot.fire(nphot_slow * visibleFraction);
            if (n > 0) {
              //calculates the time at which the photon was produced
              photon.Time = edep.T0() + GetScintTime(fRiseTimeSlow, larp->ScintSlowTimeConst(), randflatscinttime);
              // add n copies of sim::OnePhoton photon to the photon collection for a given OpChannel
              photonCollection[channel].insert(photonCollection[channel].end(), n, photon);
            }
          }
        }
      }
    }
  }
  if (lgp->UseLitePhotons()) {
    // put the photon collection of LitePhotons into the art event
    e.put(move(photLiteCol));
  }
  else {
    //put the photon collection of SimPhotons into the art event
    e.put(move(photCol));
  }
}

} // namespace phot

DEFINE_ART_MODULE(phot::PhotonLibraryPropagation)
