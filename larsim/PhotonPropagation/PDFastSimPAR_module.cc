////////////////////////////////////////////////////////////////////////
// Class:       PDFastSimPAR
// Plugin Type: producer
// File:        PDFastSimPAR_module.cc
// Description:
// - acts on sim::SimEnergyDeposit from LArG4Main,
// - simulate (fast, photon visibility service) the OpDet response to optical
// photons Input: 'sim::SimEnergyDeposit' Output: 'sim::OpDetBacktrackerRecord'
// Fast simulation of propagating the photons created from SimEnergyDeposits.

// This module does a fast simulation of propagating the photons created from
// SimEnergyDeposits, This simulation is done using the Semi-Analytical model, which
// stores the visibilities of each optical channel with respect to each optical
// voxel in the TPC volume, to avoid propagating single photons using Geant4. At
// the end of this module a collection of the propagated photons either as
//'sim::OpDetBacktrackerRecord' are placed into the art event.

// The steps this module takes are:
//  - to take number of photon and the vertex information from
//  'sim::SimEnergyDeposits',
//  - use the Semi-Analytical model (visibilities) to determine the amount of visible
//  photons at each optical channel,
//  - visible photons: the number of photons times the visibility at the middle
//  of the Geant4 step for a given optical channel.
//  - other photon information is got from 'sim::SimEnergyDeposits'
//  - add 'sim::OpDetBacktrackerRecord' to event
// Aug. 19 by Mu Wei
// Restructured Nov. 21 by P. Green
////////////////////////////////////////////////////////////////////////

// LArSoft libraries
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/IonizationScintillation/ISTPC.h"
#include "larsim/PhotonPropagation/PropagationTimeModel.h"
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"

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
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Random numbers
#include "CLHEP/Random/RandPoissonQ.h"

#include <cmath>
#include <ctime>
#include <map>
#include <memory>
#include <vector>

namespace phot {
  class PDFastSimPAR : public art::EDProducer {
  public:
    // Define the fhicl configuration
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      using DP = fhicl::DelegatedParameter;
      using ODP = fhicl::OptionalDelegatedParameter;

      fhicl::Atom<art::InputTag> SimulationLabel{Name("SimulationLabel"),
                                                 Comment("SimEnergyDeposit label.")};
      fhicl::Atom<bool> DoFastComponent{Name("DoFastComponent"),
                                        Comment("Simulate slow scintillation light, default true"),
                                        true};
      fhicl::Atom<bool> DoSlowComponent{Name("DoSlowComponent"),
                                        Comment("Simulate slow scintillation light")};
      fhicl::Atom<bool> DoReflectedLight{Name("DoReflectedLight"),
                                         Comment("Simulate reflected visible light")};
      fhicl::Atom<bool> IncludeAnodeReflections{
        Name("IncludeAnodeReflections"),
        Comment("Simulate anode reflections, default false"),
        false};
      fhicl::Atom<bool> IncludePropTime{Name("IncludePropTime"),
                                        Comment("Simulate light propagation time")};
      fhicl::Atom<bool> GeoPropTimeOnly{
        Name("GeoPropTimeOnly"),
        Comment("Simulate light propagation time geometric approximation, default false"),
        false};
      fhicl::Atom<bool> UseLitePhotons{
        Name("UseLitePhotons"),
        Comment("Store SimPhotonsLite/OpDetBTRs instead of SimPhotons")};
      fhicl::Atom<bool> OpaqueCathode{Name("OpaqueCathode"),
                                      Comment("Photons cannot cross the cathode")};
      fhicl::Atom<bool> OnlyActiveVolume{
        Name("OnlyActiveVolume"),
        Comment("PAR fast sim usually only for active volume, default true"),
        true};
      fhicl::Atom<bool> OnlyOneCryostat{Name("OnlyOneCryostat"),
                                        Comment("Set to true if light is only supported in C:1")};
      DP ScintTimeTool{Name("ScintTimeTool"),
                       Comment("Tool describing scintillation time structure")};
      ODP VUVTiming{Name("VUVTiming"), Comment("Configuration for UV timing parameterization")};
      ODP VISTiming{Name("VISTiming"),
                    Comment("Configuration for visible timing parameterization")};
      DP VUVHits{Name("VUVHits"), Comment("Configuration for UV visibility parameterization")};
      ODP VISHits{Name("VISHits"),
                  Comment("Configuration for visibile visibility parameterization")};
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit PDFastSimPAR(Parameters const& config);
    void produce(art::Event&) override;

  private:
    void Initialization();

    void detectedNumPhotons(std::vector<int>& DetectedNumPhotons,
                            const std::vector<double>& OpDetVisibilities,
                            const int NumPhotons) const;

    void AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
                     std::vector<int>& ChannelMap,
                     const sim::OpDetBacktrackerRecord& btr) const;

    bool isOpDetInSameTPC(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) const;

    // ISTPC
    larg4::ISTPC fISTPC;

    // semi-analytical model
    std::unique_ptr<SemiAnalyticalModel> fVisibilityModel;

    // propagation time model
    std::unique_ptr<PropagationTimeModel> fPropTimeModel;

    // random numbers
    CLHEP::HepRandomEngine& fPhotonEngine;
    std::unique_ptr<CLHEP::RandPoissonQ> fRandPoissPhot;
    CLHEP::HepRandomEngine& fScintTimeEngine;

    size_t fNOpChannels; // Pulled from geom during Initialization()

    // geometry properties
    std::vector<geo::BoxBoundedGeo> fActiveVolumes;
    int fNTPC;

    // optical detector properties
    std::vector<geo::Point_t> fOpDetCenter;

    //////////////////////
    // Input Parameters //
    //////////////////////

    // Module behavior
    const art::InputTag simTag;
    const bool fDoFastComponent;
    const bool fDoSlowComponent;
    const bool fDoReflectedLight;
    const bool fIncludeAnodeReflections;
    const bool fIncludePropTime;
    const bool fGeoPropTimeOnly;
    const bool fUseLitePhotons;
    const bool fOpaqueCathode;
    const bool fOnlyActiveVolume;
    const bool fOnlyOneCryostat;
    std::unique_ptr<ScintTime> fScintTime; // Tool to retrive timinig of scintillation

    // Parameterized Simulation
    fhicl::ParameterSet fVUVTimingParams;
    fhicl::ParameterSet fVISTimingParams;
    fhicl::ParameterSet fVUVHitsParams;
    fhicl::ParameterSet fVISHitsParams;
  };

  //......................................................................
  PDFastSimPAR::PDFastSimPAR(Parameters const& config)
    : art::EDProducer{config}
    , fISTPC{*(lar::providerFrom<geo::Geometry>())}
    , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,
                                                                              "HepJamesRandom",
                                                                              "photon",
                                                                              config.get_PSet(),
                                                                              "SeedPhoton"))
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,
                                                                                 "HepJamesRandom",
                                                                                 "scinttime",
                                                                                 config.get_PSet(),
                                                                                 "SeedScintTime"))
    , simTag(config().SimulationLabel())
    , fDoFastComponent(config().DoFastComponent())
    , fDoSlowComponent(config().DoSlowComponent())
    , fDoReflectedLight(config().DoReflectedLight())
    , fIncludeAnodeReflections(config().IncludeAnodeReflections())
    , fIncludePropTime(config().IncludePropTime())
    , fGeoPropTimeOnly(config().GeoPropTimeOnly())
    , fUseLitePhotons(config().UseLitePhotons())
    , fOpaqueCathode(config().OpaqueCathode())
    , fOnlyActiveVolume(config().OnlyActiveVolume())
    , fOnlyOneCryostat(config().OnlyOneCryostat())
    , fScintTime{art::make_tool<phot::ScintTime>(config().ScintTimeTool.get<fhicl::ParameterSet>())}
    , fVUVHitsParams(config().VUVHits.get<fhicl::ParameterSet>())
  {
    // Validate configuration options
    if (fIncludePropTime &&
        !config().VUVTiming.get_if_present<fhicl::ParameterSet>(fVUVTimingParams)) {
      throw art::Exception(art::errors::Configuration)
        << "Propagation time simulation requested, but VUVTiming not specified."
        << "\n";
    }

    if (fDoReflectedLight &&
        !config().VISHits.get_if_present<fhicl::ParameterSet>(fVISHitsParams)) {
      throw art::Exception(art::errors::Configuration)
        << "Reflected light simulation requested, but VisHits not specified."
        << "\n";
    }

    if (fDoReflectedLight && fIncludePropTime &&
        !config().VISTiming.get_if_present<fhicl::ParameterSet>(fVISTimingParams)) {
      throw art::Exception(art::errors::Configuration)
        << "Reflected light propagation time simulation requested, but VISTiming not specified."
        << "\n";
    }

    if (fIncludeAnodeReflections &&
        !config().VISHits.get_if_present<fhicl::ParameterSet>(fVISHitsParams)) {
      throw art::Exception(art::errors::Configuration)
        << "Anode reflections light simulation requested, but VisHits not specified."
        << "\n";
    }

    Initialization();

    if (fUseLitePhotons) {
      mf::LogInfo("PDFastSimPAR") << "Using Lite Photons";
      produces<std::vector<sim::SimPhotonsLite>>();
      produces<std::vector<sim::OpDetBacktrackerRecord>>();
      if (fDoReflectedLight) {
        mf::LogInfo("PDFastSimPAR") << "Storing Reflected Photons";
        produces<std::vector<sim::SimPhotonsLite>>("Reflected");
        produces<std::vector<sim::OpDetBacktrackerRecord>>("Reflected");
      }
    }
    else {
      mf::LogInfo("PDFastSimPAR") << "Using Sim Photons";
      produces<std::vector<sim::SimPhotons>>();
      if (fDoReflectedLight) {
        mf::LogInfo("PDFastSimPAR") << "Storing Reflected Photons";
        produces<std::vector<sim::SimPhotons>>("Reflected");
      }
    }
  }

  //......................................................................
  void PDFastSimPAR::produce(art::Event& event)
  {
    mf::LogTrace("PDFastSimPAR") << "PDFastSimPAR Module Producer"
                                 << "EventID: " << event.event();

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
      mf::LogError("PDFastSimPAR") << "PDFastSimPAR Module Cannot getByLabel: " << simTag;
      return;
    }

    auto const& edeps = edepHandle;

    int num_points = 0;
    int num_fastph = 0;
    int num_slowph = 0;
    int num_fastdp = 0;
    int num_slowdp = 0;

    for (auto const& edepi : *edeps) {
      num_points++;

      int nphot_fast = edepi.NumFPhotons();
      int nphot_slow = edepi.NumSPhotons();

      num_fastph += nphot_fast;
      num_slowph += nphot_slow;

      if (!((nphot_fast > 0 && fDoFastComponent) || (nphot_slow > 0 && fDoSlowComponent))) continue;

      int trackID = edepi.TrackID();
      int nphot = edepi.NumPhotons();
      double edeposit = edepi.Energy() / nphot;
      double pos[3] = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
      geo::Point_t const ScintPoint = {pos[0], pos[1], pos[2]};

      if (fOnlyActiveVolume && !fISTPC.isScintInActiveVolume(ScintPoint)) continue;

      // direct light
      std::vector<int> DetectedNumFast(fNOpChannels);
      std::vector<int> DetectedNumSlow(fNOpChannels);

      std::vector<double> OpDetVisibilities;
      fVisibilityModel->detectedDirectVisibilities(OpDetVisibilities, ScintPoint);
      detectedNumPhotons(DetectedNumFast, OpDetVisibilities, nphot_fast);
      detectedNumPhotons(DetectedNumSlow, OpDetVisibilities, nphot_slow);

      if (fIncludeAnodeReflections) {
        std::vector<int> AnodeDetectedNumFast(fNOpChannels);
        std::vector<int> AnodeDetectedNumSlow(fNOpChannels);

        std::vector<double> OpDetVisibilitiesAnode;
        fVisibilityModel->detectedReflectedVisibilities(OpDetVisibilitiesAnode, ScintPoint, true);
        detectedNumPhotons(AnodeDetectedNumFast, OpDetVisibilitiesAnode, nphot_fast);
        detectedNumPhotons(AnodeDetectedNumSlow, OpDetVisibilitiesAnode, nphot_slow);

        // add to existing count
        for (size_t i = 0; i < AnodeDetectedNumFast.size(); ++i) {
          DetectedNumFast[i] += AnodeDetectedNumFast[i];
        }
        for (size_t i = 0; i < AnodeDetectedNumSlow.size(); ++i) {
          DetectedNumSlow[i] += AnodeDetectedNumSlow[i];
        }
      }

      // reflected light, if enabled
      std::vector<int> ReflDetectedNumFast(fNOpChannels);
      std::vector<int> ReflDetectedNumSlow(fNOpChannels);
      if (fDoReflectedLight) {
        std::vector<double> OpDetVisibilitiesRefl;
        fVisibilityModel->detectedReflectedVisibilities(OpDetVisibilitiesRefl, ScintPoint, false);
        detectedNumPhotons(ReflDetectedNumFast, OpDetVisibilitiesRefl, nphot_fast);
        detectedNumPhotons(ReflDetectedNumSlow, OpDetVisibilitiesRefl, nphot_slow);
      }

      // loop through direct photons then reflected photons cases
      size_t DoReflected = (fDoReflectedLight) ? 1 : 0;
      for (size_t Reflected = 0; Reflected <= DoReflected; ++Reflected) {
        for (size_t channel = 0; channel < fNOpChannels; channel++) {

          if (fOpaqueCathode && !isOpDetInSameTPC(ScintPoint, fOpDetCenter[channel])) continue;

          int ndetected_fast = DetectedNumFast[channel];
          int ndetected_slow = DetectedNumSlow[channel];
          if (Reflected) {
            ndetected_fast = ReflDetectedNumFast[channel];
            ndetected_slow = ReflDetectedNumSlow[channel];
          }
          if (!((ndetected_fast > 0 && fDoFastComponent) ||
                (ndetected_slow > 0 && fDoSlowComponent)))
            continue;

          // calculate propagation time, does not matter whether fast or slow photon
          std::vector<double> transport_time;
          if (fIncludePropTime) {
            transport_time.resize(ndetected_fast + ndetected_slow);
            fPropTimeModel->propagationTime(transport_time, ScintPoint, channel, Reflected);
          }

          // SimPhotonsLite case
          if (fUseLitePhotons) {
            sim::OpDetBacktrackerRecord tmpbtr(channel);
            if (ndetected_fast > 0 && fDoFastComponent) {
              int n = ndetected_fast;
              num_fastdp += n;
              for (int i = 0; i < n; ++i) {
                // calculates the time at which the photon was produced
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
            if (ndetected_slow > 0 && fDoSlowComponent) {
              int n = ndetected_slow;
              num_slowdp += n;
              for (int i = 0; i < n; ++i) {
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
              int n = ndetected_fast;
              num_fastdp += n;
              for (int i = 0; i < n; ++i) {
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
              int n = ndetected_slow;
              num_slowdp += n;
              for (int i = 0; i < n; ++i) {
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

    mf::LogTrace("PDFastSimPAR") << "Total points: " << num_points
                                 << ", total fast photons: " << num_fastph
                                 << ", total slow photons: " << num_slowph
                                 << "\ndetected fast photons: " << num_fastdp
                                 << ", detected slow photons: " << num_slowdp;

    if (fUseLitePhotons) {
      event.put(move(phlit));
      event.put(move(opbtr));
      if (fDoReflectedLight) {
        event.put(move(phlit_ref), "Reflected");
        event.put(move(opbtr_ref), "Reflected");
      }
    }
    else {
      event.put(move(phot));
      if (fDoReflectedLight) { event.put(move(phot_ref), "Reflected"); }
    }

    return;
  }

  //......................................................................
  void PDFastSimPAR::AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
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

  //......................................................................
  void PDFastSimPAR::Initialization()
  {
    std::cout << "PDFastSimPAR Initialization" << std::endl;
    std::cout << "Initializing the geometry of the detector." << std::endl;
    std::cout << "Simulate using semi-analytic model for number of hits." << std::endl;

    fRandPoissPhot = std::make_unique<CLHEP::RandPoissonQ>(fPhotonEngine);
    geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());

    // photo-detector visibility model (semi-analytical model)
    fVisibilityModel = std::make_unique<SemiAnalyticalModel>(
      fVUVHitsParams, fVISHitsParams, fDoReflectedLight, fIncludeAnodeReflections);

    // propagation time model
    if (fIncludePropTime)
      fPropTimeModel = std::make_unique<PropagationTimeModel>(
        fVUVTimingParams, fVISTimingParams, fScintTimeEngine, fDoReflectedLight, fGeoPropTimeOnly);

    // Store info from the Geometry service
    fNOpChannels = geom.NOpDets();
    fActiveVolumes = fISTPC.extractActiveLArVolume(geom);
    fNTPC = geom.NTPC();

    {
      auto log = mf::LogTrace("PDFastSimPAR") << "PDFastSimPAR: active volume boundaries from "
                                              << fActiveVolumes.size() << " volumes:";
      for (auto const& [iCryo, box] : util::enumerate(fActiveVolumes)) {
        log << "\n - C:" << iCryo << ": " << box.Min() << " -- " << box.Max() << " cm";
      }
    } // local scope

    if (geom.Ncryostats() > 1U) {
      if (fOnlyOneCryostat) {
        mf::LogWarning("PDFastSimPAR")
          << std::string(80, '=') << "\nA detector with " << geom.Ncryostats()
          << " cryostats is configured"
          << " , and semi-analytic model is requested for scintillation photon propagation."
          << " THIS CONFIGURATION IS NOT SUPPORTED and it is open to bugs"
          << " (e.g. scintillation may be detected only in cryostat #0)."
          << "\nThis would be normally a fatal error, but it has been forcibly overridden."
          << "\n"
          << std::string(80, '=');
      }
      else {
        throw art::Exception(art::errors::Configuration)
          << "Photon propagation via semi-analytic model is not supported yet"
          << " on detectors with more than one cryostat.";
      }
    }

    for (size_t const i : util::counter(fNOpChannels)) {
      geo::OpDetGeo const& opDet = geom.OpDetGeoFromOpDet(i);
      fOpDetCenter.push_back(opDet.GetCenter());
    }
  }

  //......................................................................
  // calculates number of photons detected given visibility and emitted number of photons
  void PDFastSimPAR::detectedNumPhotons(std::vector<int>& DetectedNumPhotons,
                                        const std::vector<double>& OpDetVisibilities,
                                        const int NumPhotons) const
  {
    for (size_t i = 0; i < OpDetVisibilities.size(); ++i) {
      DetectedNumPhotons[i] = fRandPoissPhot->fire(OpDetVisibilities[i] * NumPhotons);
    }
  }

  //......................................................................
  // checks whether photo-detector is able to see the emitted light scintillation
  bool PDFastSimPAR::isOpDetInSameTPC(geo::Point_t const& ScintPoint,
                                      geo::Point_t const& OpDetPoint) const
  {
    // check optical channel is in same TPC as scintillation light, if not doesn't see light
    // temporary method working for SBND, uBooNE, DUNE 1x2x6; to be replaced to work in full DUNE geometry
    // check x coordinate has same sign or is close to zero
    if (((ScintPoint.X() < 0.) != (OpDetPoint.X() < 0.)) && std::abs(OpDetPoint.X()) > 10. &&
        fNTPC == 2) { // TODO: replace with geometry service method
      return false;
    }
    return true;
  }

  // ---------------------------------------------------------------------------

} // namespace phot

DEFINE_ART_MODULE(phot::PDFastSimPAR)
