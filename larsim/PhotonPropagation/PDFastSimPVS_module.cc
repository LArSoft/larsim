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

// LArSoft libraries
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/PropagationTimeModel.h"
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "nurandom/RandomUtils/NuRandomService.h"

// Art libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace phot {
  class PDFastSimPVS : public art::EDProducer {
  public:
    explicit PDFastSimPVS(fhicl::ParameterSet const&);
    void produce(art::Event&) override;
    void AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
                     std::vector<int>& ChannelMap,
                     const sim::OpDetBacktrackerRecord& btr) const;

  private:
    art::ServiceHandle<PhotonVisibilityService const> fPVS;
    const bool fDoFastComponent;
    const bool fDoSlowComponent;
    const bool fIncludePropTime;
    const bool fUseLitePhotons;
    const bool fStoreReflected;
    const size_t fNOpChannels;
    const art::InputTag simTag;
    CLHEP::HepRandomEngine& fScintTimeEngine;
    std::unique_ptr<ScintTime> fScintTime; // Tool to retrive timing of scintillation
    CLHEP::HepRandomEngine& fPhotonEngine;
    std::unique_ptr<CLHEP::RandPoissonQ> fRandPoissPhot;

    std::unique_ptr<PropagationTimeModel> fPropTimeModel;

    fhicl::ParameterSet fVUVTimingParams;
    fhicl::ParameterSet fVISTimingParams;
  };

  //......................................................................
  PDFastSimPVS::PDFastSimPVS(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fPVS(art::ServiceHandle<PhotonVisibilityService const>())
    , fDoFastComponent(pset.get<bool>("DoFastComponent", true))
    , fDoSlowComponent(pset.get<bool>("DoSlowComponent", true))
    , fIncludePropTime(pset.get<bool>("IncludePropTime", false))
    , fUseLitePhotons(art::ServiceHandle<sim::LArG4Parameters const>()->UseLitePhotons())
    , fStoreReflected(fPVS->StoreReflected())
    , fNOpChannels(fPVS->NOpChannels())
    , simTag{pset.get<art::InputTag>("SimulationLabel")}
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,
                                                                                 "HepJamesRandom",
                                                                                 "scinttime",
                                                                                 pset,
                                                                                 "SeedScintTime"))
    , fScintTime{art::make_tool<phot::ScintTime>(pset.get<fhicl::ParameterSet>("ScintTimeTool"))}
    , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService> {}
                      ->createEngine(*this, "HepJamesRandom", "photon", pset, "SeedPhoton"))
    , fRandPoissPhot(std::make_unique<CLHEP::RandPoissonQ>(fPhotonEngine))
  {
    mf::LogInfo("PDFastSimPVS") << "PDFastSimPVS Module Construct";

    if (fUseLitePhotons) {
      mf::LogInfo("PDFastSimPVS") << "Use Lite Photon." << std::endl;
      produces<std::vector<sim::SimPhotonsLite>>();
      produces<std::vector<sim::OpDetBacktrackerRecord>>();

      if (fStoreReflected) {
        mf::LogInfo("PDFastSimPVS") << "Store Reflected Photons";
        produces<std::vector<sim::SimPhotonsLite>>("Reflected");
        produces<std::vector<sim::OpDetBacktrackerRecord>>("Reflected");
      }
    }
    else {
      mf::LogInfo("PDFastSimPVS") << "Use Sim Photon.";
      produces<std::vector<sim::SimPhotons>>();
      if (fStoreReflected) {
        mf::LogInfo("PDFastSimPVS") << "Store Reflected Photons";
        produces<std::vector<sim::SimPhotons>>("Reflected");
      }
    }

    // Propagation times
    // validate configuration
    if (fIncludePropTime &&
        !pset.get_if_present<fhicl::ParameterSet>("VUVTiming", fVUVTimingParams)) {
      throw art::Exception(art::errors::Configuration)
        << "Propagation time simulation requested, but VUVTiming not specified."
        << "\n";
    }

    if (fIncludePropTime && fStoreReflected &&
        !pset.get_if_present<fhicl::ParameterSet>("VISTiming", fVISTimingParams)) {
      throw art::Exception(art::errors::Configuration)
        << "Reflected light propagation time simulation requested, but VISTiming not specified."
        << "\n";
    }

    // construct propagation time class
    if (fIncludePropTime)
      fPropTimeModel = std::make_unique<PropagationTimeModel>(
        fVUVTimingParams, fVISTimingParams, fScintTimeEngine, fStoreReflected);
  }

  //......................................................................
  void PDFastSimPVS::produce(art::Event& event)
  {
    std::vector<int> PDChannelToSOCMapDirect(fNOpChannels, -1);  // Where each OpChan is.
    std::vector<int> PDChannelToSOCMapReflect(fNOpChannels, -1); // Where each OpChan is.

    // SimPhotonsLite
    auto phlit = std::make_unique<std::vector<sim::SimPhotonsLite>>();
    auto opbtr = std::make_unique<std::vector<sim::OpDetBacktrackerRecord>>();
    auto phlit_ref = std::make_unique<std::vector<sim::SimPhotonsLite>>();
    auto opbtr_ref = std::make_unique<std::vector<sim::OpDetBacktrackerRecord>>();
    auto& dir_phlitcol(*phlit);
    auto& ref_phlitcol(*phlit_ref);
    // SimPhotons
    auto phot = std::make_unique<std::vector<sim::SimPhotons>>();
    auto phot_ref = std::make_unique<std::vector<sim::SimPhotons>>();
    auto& dir_photcol(*phot);
    auto& ref_photcol(*phot_ref);
    if (fUseLitePhotons) {
      dir_phlitcol.resize(fNOpChannels);
      ref_phlitcol.resize(fNOpChannels);
      for (unsigned int i = 0; i < fNOpChannels; i++) {
        dir_phlitcol[i].OpChannel = i;
        ref_phlitcol[i].OpChannel = i;
      }
    }
    else { // SimPhotons
      dir_photcol.resize(fNOpChannels);
      ref_photcol.resize(fNOpChannels);
      for (unsigned int i = 0; i < fNOpChannels; i++) {
        dir_photcol[i].fOpChannel = i;
        ref_photcol[i].fOpChannel = i;
      }
    }

    art::Handle<std::vector<sim::SimEnergyDeposit>> edepHandle;
    if (!event.getByLabel(simTag, edepHandle)) {
      mf::LogWarning("PDFastSimPVS") << "PDFastSimPVS Module Cannot getByLabel: " << simTag;
      return;
    }

    int num_points = 0;
    auto const& edeps = edepHandle;
    for (auto const& edepi : *edeps) {
      num_points++;

      int nphot_fast = edepi.NumFPhotons();
      int nphot_slow = edepi.NumSPhotons();
      if (!((nphot_fast > 0 && fDoFastComponent) || (nphot_slow > 0 && fDoSlowComponent))) continue;

      auto const& prt = edepi.MidPoint();

      MappedCounts_t Visibilities = fPVS->GetAllVisibilities(prt);
      MappedCounts_t Visibilities_Ref;
      if (fStoreReflected) {
        Visibilities_Ref = fPVS->GetAllVisibilities(prt, true);
        if (!Visibilities_Ref)
          mf::LogWarning("PDFastSimPVS") << "Fail to get visibilities for reflected photons.";
      }

      if (!Visibilities) {
        //throw cet::exception("PDFastSimPVS")
        mf::LogWarning("PDFastSimPVS")
          << "There is no entry in the PhotonLibrary for this position in space. Position: "
          << edepi.MidPoint() << "\n Move to next point";
        continue;
      }

      int trackID = edepi.TrackID();
      int nphot = edepi.NumPhotons();
      double edeposit = edepi.Energy() / nphot;
      double pos[3] = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};

      geo::Point_t const ScintPoint = {pos[0], pos[1], pos[2]};

      // loop through direct photons then reflected photons cases
      size_t DoReflected = (fStoreReflected) ? 1 : 0;
      for (size_t Reflected = 0; Reflected <= DoReflected; ++Reflected) {
        for (unsigned int channel = 0; channel < fNOpChannels; ++channel) {

          double visibleFraction = (Reflected) ? Visibilities_Ref[channel] : Visibilities[channel];
          if (visibleFraction < 1e-9) continue; // voxel is not visible at this optical channel

          // number of detected photons
          int ndetected_fast =
            (nphot_fast > 0) ? fRandPoissPhot->fire(nphot_fast * visibleFraction) : 0;
          int ndetected_slow =
            (nphot_slow > 0) ? fRandPoissPhot->fire(nphot_slow * visibleFraction) : 0;
          if (!((ndetected_fast > 0 && fDoFastComponent) ||
                (ndetected_slow > 0 && fDoSlowComponent)))
            continue;

          // calculate propagation times if included, does not matter whether fast or slow photon
          std::vector<double> transport_time;
          if (fIncludePropTime) {
            transport_time.resize(ndetected_fast + ndetected_slow);
            fPropTimeModel->propagationTime(transport_time, ScintPoint, channel, Reflected);
          }

          // SimPhotonsLite case
          if (fUseLitePhotons) {
            sim::OpDetBacktrackerRecord tmpbtr(channel);
            if (ndetected_fast > 0 && fDoFastComponent) {
              for (int i = 0; i < ndetected_fast; ++i) {
                // calculate the time at which each photon is seen
                fScintTime->GenScintTime(true, fScintTimeEngine);
                int time;
                if (fIncludePropTime)
                  time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() +
                                          transport_time[i]);
                else
                  time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                if (Reflected)
                  ++ref_phlitcol[channel].DetectedPhotons[time];
                else
                  ++dir_phlitcol[channel].DetectedPhotons[time];
                tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
              }
            }
            if ((ndetected_slow > 0) && fDoSlowComponent) {
              for (int i = 0; i < ndetected_slow; ++i) {
                // calculate the time at which each photon is seen
                fScintTime->GenScintTime(false, fScintTimeEngine);
                int time;
                if (fIncludePropTime)
                  time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() +
                                          transport_time[ndetected_fast + i]);
                else
                  time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                if (Reflected)
                  ++ref_phlitcol[channel].DetectedPhotons[time];
                else
                  ++dir_phlitcol[channel].DetectedPhotons[time];
                tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
              }
            }
            if (Reflected)
              AddOpDetBTR(*opbtr_ref, PDChannelToSOCMapReflect, tmpbtr);
            else
              AddOpDetBTR(*opbtr, PDChannelToSOCMapDirect, tmpbtr);
          }

          // SimPhotons case
          else {
            sim::OnePhoton photon;
            photon.SetInSD = false;
            photon.InitialPosition = edepi.End();
            if (Reflected)
              photon.Energy = 2.9 * CLHEP::eV; // 430 nm
            else
              photon.Energy = 9.7 * CLHEP::eV; // 128 nm
            if (ndetected_fast > 0 && fDoFastComponent) {
              for (int i = 0; i < ndetected_fast; ++i) {
                // calculates the time at which the photon was produced
                fScintTime->GenScintTime(true, fScintTimeEngine);
                int time;
                if (fIncludePropTime)
                  time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() +
                                          transport_time[i]);
                else
                  time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                photon.Time = time;
                if (Reflected)
                  ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
                else
                  dir_photcol[channel].insert(dir_photcol[channel].end(), 1, photon);
              }
            }
            if (ndetected_slow > 0 && fDoSlowComponent) {
              for (int i = 0; i < ndetected_slow; ++i) {
                fScintTime->GenScintTime(false, fScintTimeEngine);
                int time;
                if (fIncludePropTime)
                  time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() +
                                          transport_time[ndetected_fast + i]);
                else
                  time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                photon.Time = time;
                if (Reflected)
                  ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
                else
                  dir_photcol[channel].insert(dir_photcol[channel].end(), 1, photon);
              }
            }
          }
        }
      }
    }

    if (fUseLitePhotons) {
      event.put(move(phlit));
      event.put(move(opbtr));
      if (fStoreReflected) {
        event.put(move(phlit_ref), "Reflected");
        event.put(move(opbtr_ref), "Reflected");
      }
    }
    else {
      event.put(move(phot));
      if (fStoreReflected) { event.put(move(phot_ref), "Reflected"); }
    }

    return;
  }

  //......................................................................
  void PDFastSimPVS::AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
                                 std::vector<int>& ChannelMap,
                                 const sim::OpDetBacktrackerRecord& btr) const
  {
    int iChan = btr.OpDetNum();
    if (ChannelMap[iChan] < 0) {
      ChannelMap[iChan] = opbtr.size();
      opbtr.emplace_back(std::move(btr));
    }
    else {
      size_t idtest = ChannelMap[iChan];
      auto const& timePDclockSDPsMap = btr.timePDclockSDPsMap();
      for (auto const& timePDclockSDP : timePDclockSDPsMap) {
        for (auto const& sdp : timePDclockSDP.second) {
          double xyz[3] = {sdp.x, sdp.y, sdp.z};
          opbtr.at(idtest).AddScintillationPhotons(
            sdp.trackID, timePDclockSDP.first, sdp.numPhotons, xyz, sdp.energy);
        }
      }
    }
  }

} // namespace

DEFINE_ART_MODULE(phot::PDFastSimPVS)
