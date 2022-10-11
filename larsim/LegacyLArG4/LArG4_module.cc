////////////////////////////////////////////////////////////////////////
/// \file  LArG4_module.cc
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// This a module.  It has the following functions:
///
/// - Initialize Geant4 physics, detector geometry, and other
///   processing.
///
/// - Accept sim::MCTruth objects from the MC branch of the FMWK
///   Event structure.
///
/// - Pass the primary particles to the Geant4 simulation to calculate
///   "truth" information for the detector response.
///
/// - Pass the truth information to the DetSim branch of the FMWK event.

#include "nug4/G4Base/G4Helper.h"

// C++ Includes
#include <cassert>
#include <map>
#include <set>
#include <sstream>
#include <sys/stat.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/ParticleFilters.h" // util::PositionInVolumeFilter
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/MCDumpers/MCDumpers.h" // sim::dump namespace
#include "lardataobj/MCBase/MCParticleLite.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/LegacyLArG4/AllPhysicsLists.h"
#include "larsim/LegacyLArG4/AuxDetReadout.h"
#include "larsim/LegacyLArG4/AuxDetReadoutGeometry.h"
#include "larsim/LegacyLArG4/IonizationAndScintillation.h"
#include "larsim/LegacyLArG4/LArStackingAction.h"
#include "larsim/LegacyLArG4/LArVoxelReadout.h"
#include "larsim/LegacyLArG4/LArVoxelReadoutGeometry.h"
#include "larsim/LegacyLArG4/MaterialPropertyLoader.h"
#include "larsim/LegacyLArG4/OpDetPhotonTable.h"
#include "larsim/LegacyLArG4/OpDetReadoutGeometry.h"
#include "larsim/LegacyLArG4/OpDetSensitiveDetector.h"
#include "larsim/LegacyLArG4/ParticleListAction.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nug4/G4Base/UserActionManager.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// G4 Includes
#include "Geant4/G4LogicalVolumeStore.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/G4VUserDetectorConstruction.hh"

//For energy depositions
#include "lardataobj/Simulation/SimEnergyDeposit.h"

// Boost includes
#include "boost/algorithm/string.hpp"

///Geant4 interface
namespace larg4 {

  // Forward declarations within namespace.
  class LArVoxelListAction;

  /**
   * @brief Runs Geant4 simulation and propagation of electrons and photons to readout
   *
   * This module collects generated particles from one or more generators and
   * processes them through Geant4.
   *
   * Input
   * ------
   *
   * The module reads the particles to process from `simb::MCTruth` records.
   * Each particle generator is required to produce a vector of such records:
   * `std::vector<simb::MCTruth>`.
   *
   * The module allows two operation modes:
   * -# process specific generators: the label of the generator modules to be
   *   processed is specified explicitly in `LArG4` configuration
   * -# process all truth information generated so far: no generator is specified
   *   in the `LArG4` module configuration, and the module will process all
   *   data products of type `std::vector<simb::MCTruth>`, in a non-specified
   *   order
   *
   * For each `simb::MCTruth`, a Geant4 run is started.
   * The interface with Geant4 is via a helper class provided by _nug4_.
   * Only the particles in the truth record which have status code
   * (`simb::MCParticle::StatusCode()`) equal to `1` are processed.
   * These particles are called, in `LArG4` jargon, _primaries_.
   *
   *
   * Output
   * -------
   *
   * The `LArG4` module produces:
   * * a collection of `sim::SimChannel`: each `sim::SimChannel` represents the
   *   set of energy depositions in liquid argon which drifted and were observed
   *   on a certain channel; it includes physics effects like attenuation,
   *   diffusion, electric field distortion, etc. Information of the generating
   *   Geant4 "track" is retained;
   * * a collection of `sim::SimPhotons` or `sim::SimPhotonsLite`: each
   *   `sim::SimPhotons` represents the set of individual photons reaching a
   *   channel of the optical detector; it includes physics effects as well as
   *   quantum efficiency of the detector (to reduce data size early in the
   *   process); `sim::SimPhotonsLite` drops the information of the single
   *   photons and stores only collective information (e.g. their number).
   * * a collection of `sim::OpDetBacktrackerRecord` (to be documented)
   * * a collection of `sim::AuxDetSimChannel` (to be documented)
   * * a collection of `simb::MCParticle`: the particles generated in the
   *   interaction of the primary particles with the material in the world
   *   are stored, but minor filtering by geometry and by physics is possible.
   *   An association of them with the originating `simb::MCTruth` object is
   *   also produced.
   *
   *
   * Notes on the conventions
   * -------------------------
   *
   * * all and the particles in the truth record (`simb::MCTruth`) which have
   *   status code (`simb::MCParticle::StatusCode()`) equal to `1` are passed
   *   to Geant4. These particles are called, in `LArG4` jargon, _primaries_.
   *   The interface with Geant4 is via a helper class provided by _nug4_.
   * * normally, information about each particle that Geant4 propagates (which
   *   Geant4 calls _tracks_), primary or not, is saved as an individual
   *   `simb::MCParticle` object into the output particle list. Each
   *   `simb::MCParticle` includes a Geant4-like track ID which is also recorded
   *   into each `sim::IDE` deposited by that particle. This information can be
   *   used to track all the deposition from a particle, or to backtrack the
   *   particle responsible of a deposition (but see below...).
   *   Note that the stored track ID may be different than the one Geant4 used
   *   (and, in particular, it's guaranteed to be unique within a `sim::LArG4`
   *   instance output).
   * * there are options (some set in `sim::LArG4Parameters` service) which
   *   allow for Geant4 tracks not to be saved as `simb::MCParticle` (e.g.
   *   `ParticleKineticEnergyCut`, `KeepEMShowerDaughters`). When these
   *   particles have deposited energy, their `sim::IDE` will report the ID of
   *   the first parent Geant4 track which is saved in the `simb::MCParticle`
   *   list, but _with its sign flipped_. Therefore, when tracking or
   *   backtracking (see above), comparisons should be performed using the
   *   absolute value of the `sim::IDE` (e.g. `std::abs(ide.trackID)`).
   *
   *
   * Timing
   * -------
   *
   * The `LArG4` module produces `sim::SimChannel` objects from generated
   * `simb::MCParticle`. Each particle ("primary") is assigned the time taken
   * from its vertex (a 4-vector), which is expected to be represented in
   * nanoseconds.
   * The `sim::SimChannel` object is a collection of `sim::IDE` in time. The
   * position in the `sim::IDE` is the location where some ionization occurred.
   * The time associated to a `sim::IDE` is stored in tick units. The time it
   * represents is the time when the ionization happened, which is the time of
   * the primary particle plus the propagation time to the ionization location,
   * plus the drift time, which the ionized electrons take to reach the anode
   * wire. This time is then shifted to the frame of the electronics time
   * via `detinfo::DetectorClocks::G4ToElecTime()`, which adds a configurable
   * time offset. The time is converted into ticks via
   * `detinfo::DetectorClocks::TPCClock()`, and this is the final value
   * associated to the `sim::IDE`. For a more complete overview, see
   * https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/Simulation#Simulation-Timing
   *
   *
   * Randomness
   * -----------
   *
   * The random number generators used by this process are:
   * - 'GEANT' instance: used by Geant4
   * - 'propagation' instance: used in electron propagation
   *
   *
   * Configuration parameters
   * -------------------------
   *
   * - *G4PhysListName* (string, default: `"larg4::PhysicsList"`):
   *     whether to use the G4 overlap checker, which catches different issues than ROOT
   * - *CheckOverlaps* (bool, default: `false`):
   *     whether to use the G4 overlap checker
   * - *DumpParticleList* (bool, default: `false`):
   *     whether to print all MCParticles tracked;
   *     requires `MakeMCParticles` being `true`
   * - *DumpSimChannels* (bool, default: `false`):
   *     whether to print all depositions on each SimChannel
   * - *SmartStacking* (int, default: `0`):
   *     whether to use class to dictate how tracks are put on stack (nonzero is on)
   * - *MakeMCParticles* (flag, default: `true`): keep a list of the particles
   *     seen in the detector, and eventually save it; you almost always want this on
   * - *KeepParticlesInVolumes* (list of strings, default: _empty_):
   *     list of volumes in which to keep `simb::MCParticle` objects (empty keeps all);
   *     requires `MakeMCParticles` being `true`
   * - *GeantCommandFile* (string, _required_):
   *     G4 macro file to pass to `G4Helper` for setting G4 command
   * - *Seed* (integer, not defined by default): if defined, override the seed for
   *     random number generator used in Geant4 simulation (which is obtained from
   *     `NuRandomService` by default)
   * - *PropagationSeed* (integer, not defined by default): if defined,
   *     override the seed for the random generator used for electrons propagation
   *     to the wire planes (obtained from the `NuRandomService` by default)
   * - *InputLabels* (list of strings, default: process all truth):
   *     optional list of generator labels whose produced `simb::MCTruth` will
   *     be simulated; if not specified, all `simb::MCTruth` vector data
   *     products are simulated
   * - *ChargeRecoveryMargin* (double, default: `0`): sets the maximum
   *     distance from a plane for the wire charge recovery to occur, in
   *     centimeters; for details on how it works, see
   *     `larg4::LArVoxelReadout::SetOffPlaneChargeRecoveryMargin()`. A value of
   *     `0` effectively disables this feature. All TPCs will have the same
   *     margin applied.
   *
   *
   * Simulation details
   * ===================
   *
   * Source of the operational parameters
   * -------------------------------------
   *
   * @anchor LArG4_MaterialProperties
   *
   * Some of the physical properties have their values set in FHiCL
   * configuration (e.g. `detinfo::LArParameters`). Then, GEANT4 is informed
   * of them via `larg4::MaterialPropertyLoader`. The material property table
   * in GEANT4 is then used by other LArSoft components to discover the
   * parameter values.
   *
   * Among the parameters registered to GEANT4, the scintillation yields, i.e.
   * how many scintillation photons are produced on average by 1 MeV of
   * deposited energy, are also stored by type of ioniziong particle.
   * These scintillation yields _do include a prescale factor_ (that may
   * include, for example, the photomultiplier quantum efficiency), from the
   * `ScintPreScale` parameter of `detinfo::LArPropertiesStandard`
   * or equivalent.
   *
   *
   * Reflectivity to optical photons
   * --------------------------------
   *
   * Two models are supported for the simulation of (scintillation) light
   * crossing detector surfaces:
   * -# the standard one from GEANT4, implemented in `G4OpBoundaryProcess`
   * -# a simplified one, implemented in `larg4::OpBoundaryProcessSimple`
   *
   * The model is chosen according to the value of
   * `detinfo::DetectorProperties::SimpleBoundary()`, and the choice is
   * currently exerted by `larg4::OpticalPhysics`.
   *
   * The simplified model is faster and simpler: it only deals with absorption
   * and reflection (both specular and diffues).
   * This is the "default" model used in most contexts.
   *
   * GEANT4 model is more complete and slower. It may take some art to fully
   * configure all the properties of the materials at the sides of the surfaces.
   * The price is a detailed simulation that includes among others refraction
   * and wavelength shifting.
   *
   *
   * Scintillation
   * --------------
   *
   * When using the fast optical simulation, which is the "standard" running
   * mode, energy depositions from GEANT4 are "converted" into a number of
   * scintillation photons by the global `larg4::IonizationAndScintillation`
   * object instance, which internally utilizes the algorithm set up via
   * configuration parameter `IonAndScintCalculator` in `LArG4Parameters`
   * service (at the time of writing, `"Separate"` is supported and `"NEST"` is
   * accepted too).
   * The number of scintillation photons per energy unit is read from GEANT4
   * @ref LArG4_MaterialProperties "material properties table". It includes
   * already quantum efficiency ("prescale") and it may depend on the type of
   * ionizing particle, depending on the configuration (`LArPropertiesStandard`
   * parameter `ScintByParticleType`). This value ("yield") is used as
   * the average of a Poisson distribution from which the actual number of
   * scintillation photons is extracted case by case.
   * The implementation `larg4::ISCalculationSeparate` may also include medium
   * saturation effects as well, if configured, but only if the scintillation
   * yield is set not to depend on the type of ionizing particle.
   * The number of scintillation photons is then distributed between the fast
   * and slow component by a yield ratio also set in the material parameters,
   * and the single photons are distributed in time accordingly to their
   * component.
   *
   */
  class LArG4 : public art::EDProducer {
  public:
    explicit LArG4(fhicl::ParameterSet const& pset);

  private:
    /// The main routine of this module: Fetch the primary particles
    /// from the event, simulate their evolution in the detector, and
    /// produce the detector response.
    void produce(art::Event& evt) override;
    void beginJob() override;
    void beginRun(art::Run& run) override;

    std::unique_ptr<g4b::G4Helper> fG4Help{nullptr}; ///< G4 interface object
    larg4::ParticleListAction* fparticleListAction{
      nullptr}; ///< Geant4 user action to particle information.

    std::string fG4PhysListName; ///< predefined physics list to use if not making a custom one
    std::string fG4MacroPath;    ///< directory path for Geant4 macro file to be
                                 ///< executed before main MC processing.
    bool fCheckOverlaps;         ///< Whether to use the G4 overlap checker
    bool fMakeMCParticles;       ///< Whether to keep a `sim::MCParticle` list
    bool
      fStoreDroppedMCParticles; ///< Whether to keep a `sim::MCParticleLite` list of dropped particles
    bool fdumpParticleList;     ///< Whether each event's sim::ParticleList will be displayed.
    bool fdumpSimChannels;      ///< Whether each event's sim::Channel will be displayed.
    bool fUseLitePhotons;
    bool fStoreReflected{false};
    int fSmartStacking;          ///< Whether to instantiate and use class to
    double fOffPlaneMargin = 0.; ///< Off-plane charge recovery margin
                                 ///< dictate how tracks are put on stack.
    std::vector<std::string> fInputLabels;
    std::vector<std::string>
      fKeepParticlesInVolumes; ///<Only write particles that have trajectories through these volumes

    bool fSparsifyTrajectories; ///< Sparsify MCParticle Trajectories

    CLHEP::HepRandomEngine& fEngine; ///< Random-number engine for IonizationAndScintillation
                                     ///< initialization

    detinfo::DetectorPropertiesData fDetProp; ///< Must outlive fAllPhysicsLists!
    AllPhysicsLists fAllPhysicsLists;
    LArVoxelReadoutGeometry* fVoxelReadoutGeometry{
      nullptr}; /// Pointer used for correctly updating the clock data state.

    /// Configures and returns a particle filter
    std::unique_ptr<util::PositionInVolumeFilter> CreateParticleVolumeFilter(
      std::set<std::string> const& vol_names) const;
  };

} // namespace LArG4

namespace {

  // ---------------------------------------------------------------------------
  /**
   * @brief Moves data from `source` to the end of `dest`.
   * @tparam T type of values in the containers
   * @param dest collection to append data to
   * @param source collection to take data from
   * @return a reference to `dest`
   *
   * The data contained in `source` is moved (appended) to the end of `dest`.
   *
   * The data is moved from source element by element; as an exception, if
   * `dest` is still empty, the data is moved in a block.
   *
   * The `source` collection is always returned depleted as after being "moved
   * away" (i.e. after `auto temp = std::move(source);`): that means that the
   * only thing that can be done with it according to C++ standard is to
   * destruct it.
   *
   */
  template <typename T>
  std::vector<T>& append(std::vector<T>& dest, std::vector<T>&& source)
  {
    if (empty(dest))
      dest = std::move(source);
    else {
      dest.insert(dest.end(), std::move_iterator{begin(source)}, std::move_iterator{end(source)});
      source = std::vector<T>{}; // ensure the old memory is released
    }
    return dest;
  }
  // ---------------------------------------------------------------------------

} // local namespace

namespace larg4 {

  //----------------------------------------------------------------------
  // Constructor
  LArG4::LArG4(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fG4PhysListName(pset.get<std::string>("G4PhysListName", "larg4::PhysicsList"))
    , fCheckOverlaps(pset.get<bool>("CheckOverlaps", false))
    , fMakeMCParticles(pset.get<bool>("MakeMCParticles", true))
    , fStoreDroppedMCParticles(pset.get<bool>("StoreDroppedMCParticles", false))
    , fdumpParticleList(pset.get<bool>("DumpParticleList", false))
    , fdumpSimChannels(pset.get<bool>("DumpSimChannels", false))
    , fSmartStacking(pset.get<int>("SmartStacking", 0))
    , fOffPlaneMargin(pset.get<double>("ChargeRecoveryMargin", 0.0))
    , fKeepParticlesInVolumes(pset.get<std::vector<std::string>>("KeepParticlesInVolumes", {}))
    , fSparsifyTrajectories(pset.get<bool>("SparsifyTrajectories", false))
    , fEngine(art::ServiceHandle<rndm::NuRandomService> {}
                ->createEngine(*this, "HepJamesRandom", "propagation", pset, "PropagationSeed"))
    , fDetProp{art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob()}
    , fAllPhysicsLists{fDetProp}
  {
    MF_LOG_DEBUG("LArG4") << "Debug: LArG4()";
    art::ServiceHandle<art::RandomNumberGenerator const> rng;

    if (!fMakeMCParticles) { // configuration option consistency
      if (fdumpParticleList) {
        throw art::Exception(art::errors::Configuration)
          << "Option `DumpParticleList` can't be set if `MakeMCParticles` is unset.\n";
      }
      if (!fKeepParticlesInVolumes.empty()) {
        throw art::Exception(art::errors::Configuration)
          << "Option `KeepParticlesInVolumes` can't be set if `MakeMCParticles` is unset.\n";
      }
    } // if

    if (pset.has_key("Seed")) {
      throw art::Exception(art::errors::Configuration)
        << "The configuration of LArG4 module has the discontinued 'Seed' parameter.\n"
           "Seeds are now controlled by two parameters: 'GEANTSeed' and 'PropagationSeed'.";
    }
    // setup the random number service for Geant4, the "G4Engine" label is a
    // special tag setting up a global engine for use by Geant4/CLHEP;
    // obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed" or "GEANTSeed"
    // FIXME: THIS APPEARS TO BE A NO-OP; IS IT NEEDED?
    (void)art::ServiceHandle<rndm::NuRandomService>()->createEngine(
      *this, "G4Engine", "GEANT", pset, "GEANTSeed");

    //get a list of generators to use, otherwise, we'll end up looking for anything that's
    //made an MCTruth object
    bool useInputLabels =
      pset.get_if_present<std::vector<std::string>>("InputLabels", fInputLabels);
    if (!useInputLabels) fInputLabels.resize(0);

    art::ServiceHandle<sim::LArG4Parameters const> lgp;
    fUseLitePhotons = lgp->UseLitePhotons();

    if (!lgp->NoPhotonPropagation()) {
      try {
        art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
        fStoreReflected = pvs->StoreReflected();
      }
      catch (art::Exception const& e) {
        // If the service is not configured, then just keep the default
        // false for reflected light. If reflected photons are simulated
        // without PVS they will show up in the regular SimPhotons collection
        if (e.categoryCode() != art::errors::ServiceNotFound) throw;
      }

      if (!fUseLitePhotons) {
        produces<std::vector<sim::SimPhotons>>();
        if (fStoreReflected) { produces<std::vector<sim::SimPhotons>>("Reflected"); }
      }
      else {
        produces<std::vector<sim::SimPhotonsLite>>();
        produces<std::vector<sim::OpDetBacktrackerRecord>>();
        if (fStoreReflected) {
          produces<std::vector<sim::SimPhotonsLite>>("Reflected");
          produces<std::vector<sim::OpDetBacktrackerRecord>>("Reflected");
        }
      }
    }

    if (lgp->FillSimEnergyDeposits()) {
      produces<std::vector<sim::SimEnergyDeposit>>("TPCActive");
      produces<std::vector<sim::SimEnergyDeposit>>("Other");
    }

    if (fMakeMCParticles) {
      produces<std::vector<simb::MCParticle>>();
      produces<art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo>>();
    }
    if (fStoreDroppedMCParticles) { produces<std::vector<sim::MCParticleLite>>(); }
    if (!lgp->NoElectronPropagation()) produces<std::vector<sim::SimChannel>>();
    produces<std::vector<sim::AuxDetSimChannel>>();

    // constructor decides if initialized value is a path or an environment variable
    cet::search_path sp("FW_SEARCH_PATH");

    sp.find_file(pset.get<std::string>("GeantCommandFile"), fG4MacroPath);
    struct stat sb;
    if (fG4MacroPath.empty() || stat(fG4MacroPath.c_str(), &sb) != 0)
      // failed to resolve the file name
      throw cet::exception("NoG4Macro") << "G4 macro file " << fG4MacroPath << " not found!\n";
  }

  //----------------------------------------------------------------------
  void LArG4::beginJob()
  {
    fG4Help = std::make_unique<g4b::G4Helper>(fG4MacroPath, fG4PhysListName);

    if (fCheckOverlaps) fG4Help->SetOverlapCheck(true);

    art::ServiceHandle<geo::Geometry const> geom;
    fG4Help->ConstructDetector(geom->GDMLFile());

    // Get the logical volume store and assign material properties
    larg4::MaterialPropertyLoader mpl;
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    mpl.GetPropertiesFromServices(detProp);
    mpl.UpdateGeometry(G4LogicalVolumeStore::GetInstance());

    // Tell the detector about the parallel LAr voxel geometry.
    std::vector<G4VUserParallelWorld*> pworlds;

    // Intialize G4 physics and primary generator action
    fG4Help->InitPhysics();

    // create the ionization and scintillation calculator;
    // this is a singleton (!) so it does not make sense
    // to create it in LArVoxelReadoutGeometry
    IonizationAndScintillation::CreateInstance(detProp, fEngine);

    // make a parallel world for each TPC in the detector
    LArVoxelReadoutGeometry::Setup_t readoutGeomSetupData;
    readoutGeomSetupData.readoutSetup.offPlaneMargin = fOffPlaneMargin;
    readoutGeomSetupData.readoutSetup.propGen = &fEngine;

    fVoxelReadoutGeometry =
      new LArVoxelReadoutGeometry("LArVoxelReadoutGeometry", readoutGeomSetupData);
    pworlds.push_back(fVoxelReadoutGeometry);
    pworlds.push_back(
      new OpDetReadoutGeometry(geom->OpDetGeoName(), "OpDetReadoutGeometry", fUseLitePhotons));
    pworlds.push_back(new AuxDetReadoutGeometry("AuxDetReadoutGeometry"));

    fG4Help->SetParallelWorlds(pworlds);

    // moved up
    // Intialize G4 physics and primary generator action
    fG4Help->InitPhysics();

    // Use the UserActionManager to handle all the Geant4 user hooks.
    g4b::UserActionManager* uaManager = g4b::UserActionManager::Instance();

    // User-action class for accumulating LAr voxels.
    art::ServiceHandle<sim::LArG4Parameters const> lgp;

    // User-action class for accumulating particles and trajectories
    // produced in the detector.
    fparticleListAction = new larg4::ParticleListAction(lgp->ParticleKineticEnergyCut(),
                                                        lgp->StoreTrajectories(),
                                                        lgp->KeepEMShowerDaughters(),
                                                        fMakeMCParticles,
                                                        fStoreDroppedMCParticles);
    uaManager->AddAndAdoptAction(fparticleListAction);

    // UserActionManager is now configured so continue G4 initialization
    fG4Help->SetUserAction();

    // With an enormous detector with lots of rock ala LAr34 (nee LAr20)
    // we need to be smarter about stacking.
    if (fSmartStacking > 0) {
      G4UserStackingAction* stacking_action = new LArStackingAction(fSmartStacking);
      fG4Help->GetRunManager()->SetUserAction(stacking_action);
    }
  }

  void LArG4::beginRun(art::Run& run)
  {
    // prepare the filter object (null if no filtering)
    std::set<std::string> volnameset(fKeepParticlesInVolumes.begin(),
                                     fKeepParticlesInVolumes.end());
    fparticleListAction->ParticleFilter(CreateParticleVolumeFilter(volnameset));
  }

  std::unique_ptr<util::PositionInVolumeFilter> LArG4::CreateParticleVolumeFilter(
    std::set<std::string> const& vol_names) const
  {
    // if we don't have favourite volumes, don't even bother creating a filter
    if (empty(vol_names)) return {};

    auto const& geom = *art::ServiceHandle<geo::Geometry const>();

    std::vector<std::vector<TGeoNode const*>> node_paths = geom.FindAllVolumePaths(vol_names);

    // collection of interesting volumes
    util::PositionInVolumeFilter::AllVolumeInfo_t GeoVolumePairs;
    GeoVolumePairs.reserve(node_paths.size()); // because we are obsessed

    //for each interesting volume, follow the node path and collect
    //total rotations and translations
    for (size_t iVolume = 0; iVolume < node_paths.size(); ++iVolume) {
      std::vector<TGeoNode const*> path = node_paths[iVolume];

      auto pTransl = new TGeoTranslation(0., 0., 0.);
      auto pRot = new TGeoRotation();
      for (TGeoNode const* node : path) {
        TGeoTranslation thistranslate(*node->GetMatrix());
        TGeoRotation thisrotate(*node->GetMatrix());
        pTransl->Add(&thistranslate);
        *pRot = *pRot * thisrotate;
      }

      // for some reason, pRot and pTransl don't have tr and rot bits set
      // correctly make new translations and rotations so bits are set correctly
      auto pTransl2 = new TGeoTranslation(
        pTransl->GetTranslation()[0], pTransl->GetTranslation()[1], pTransl->GetTranslation()[2]);
      double phi = 0., theta = 0., psi = 0.;
      pRot->GetAngles(phi, theta, psi);
      auto pRot2 = new TGeoRotation();
      pRot2->SetAngles(phi, theta, psi);

      auto pTransf = new TGeoCombiTrans(*pTransl2, *pRot2);
      GeoVolumePairs.emplace_back(node_paths[iVolume].back()->GetVolume(), pTransf);
    }

    return std::make_unique<util::PositionInVolumeFilter>(std::move(GeoVolumePairs));
  } // CreateParticleVolumeFilter()

  void LArG4::produce(art::Event& evt)
  {
    MF_LOG_DEBUG("LArG4") << "produce()";
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    LArVoxelReadoutGeometry::Sentry const set_for_event{fVoxelReadoutGeometry, clockData, detProp};

    // loop over the lists and put the particles and voxels into the event as
    // collections
    auto scCol = std::make_unique<std::vector<sim::SimChannel>>();
    auto adCol = std::make_unique<std::vector<sim::AuxDetSimChannel>>();
    auto tpassn = fMakeMCParticles ?
                    std::make_unique<
                      art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo>>() :
                    nullptr;
    auto partCol = fMakeMCParticles ? std::make_unique<std::vector<simb::MCParticle>>() : nullptr;
    auto droppedPartCol =
      fStoreDroppedMCParticles ? std::make_unique<std::vector<sim::MCParticleLite>>() : nullptr;
    auto PhotonCol = std::make_unique<std::vector<sim::SimPhotons>>();
    auto PhotonColRefl = std::make_unique<std::vector<sim::SimPhotons>>();
    auto LitePhotonCol = std::make_unique<std::vector<sim::SimPhotonsLite>>();
    auto LitePhotonColRefl = std::make_unique<std::vector<sim::SimPhotonsLite>>();
    auto cOpDetBacktrackerRecordCol = std::make_unique<std::vector<sim::OpDetBacktrackerRecord>>();
    auto cOpDetBacktrackerRecordColRefl =
      std::make_unique<std::vector<sim::OpDetBacktrackerRecord>>();

    std::optional<art::PtrMaker<simb::MCParticle>> makeMCPartPtr;
    if (fMakeMCParticles) makeMCPartPtr.emplace(evt);

    // for energy deposits
    auto edepCol_TPCActive = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
    auto edepCol_Other = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

    // Fetch the lists of LAr voxels and particles.
    art::ServiceHandle<sim::LArG4Parameters const> lgp;
    art::ServiceHandle<geo::Geometry const> geom;

    // Clear the detected photon table
    OpDetPhotonTable::Instance()->ClearTable(geom->NOpDets());
    if (lgp->FillSimEnergyDeposits()) OpDetPhotonTable::Instance()->ClearEnergyDeposits();

    // reset the track ID offset as we have a new collection of interactions
    fparticleListAction->ResetTrackIDOffset();

    //look to see if there is any MCTruth information for this
    //event
    std::vector<art::Handle<std::vector<simb::MCTruth>>> mclists;
    if (empty(fInputLabels))
      //evt.getManyByType(mclists);
      mclists = evt.getMany<std::vector<simb::MCTruth>>();
    else {
      mclists.resize(fInputLabels.size());
      for (size_t i = 0; i < fInputLabels.size(); i++)
        evt.getByLabel(fInputLabels[i], mclists[i]);
    }

    unsigned int nGeneratedParticles = 0;

    // Need to process Geant4 simulation for each interaction separately.
    for (size_t mcl = 0; mcl < mclists.size(); ++mcl) {

      art::Handle<std::vector<simb::MCTruth>> mclistHandle = mclists[mcl];

      for (size_t m = 0; m < mclistHandle->size(); ++m) {
        art::Ptr<simb::MCTruth> mct(mclistHandle, m);

        MF_LOG_DEBUG("LArG4") << *(mct.get());

        // The following tells Geant4 to track the particles in this interaction.
        fG4Help->G4Run(mct);

        if (!partCol) continue;
        assert(tpassn);

        // receive the particle list
        sim::ParticleList particleList = fparticleListAction->YieldList();

        for (auto const& partPair : particleList) {
          simb::MCParticle& p = *(partPair.second);
          ++nGeneratedParticles;

          // if the particle has been marked as dropped, we don't save it
          // (as of LArSoft ~v5.6 this does not ever happen because
          // ParticleListAction has already taken care of deleting them)
          if (ParticleListAction::isDropped(&p)) continue;

          sim::GeneratedParticleInfo const truthInfo{
            fparticleListAction->GetPrimaryTruthIndex(p.TrackId())};
          if (!truthInfo.hasGeneratedParticleIndex() && (p.Mother() == 0)) {
            // this means it's primary but with no information; logic error!!
            art::Exception error(art::errors::LogicError);
            error << "Failed to match primary particle:\n";
            sim::dump::DumpMCParticle(error, p, "  ");
            error << "\nwith particles from the truth record '"
                  << mclistHandle.provenance()->inputTag() << "':\n";
            sim::dump::DumpMCTruth(error, *mct, 2U, "  "); // 2 points per line
            error << "\n";
            throw error;
          }

          if (fSparsifyTrajectories) p.SparsifyTrajectory();

          partCol->push_back(std::move(p));

          tpassn->addSingle(mct, (*makeMCPartPtr)(partCol->size() - 1), truthInfo);

        } // for(particleList)

        if (fStoreDroppedMCParticles && droppedPartCol) {
          // Request a list of dropped particles
          // Store them in MCParticleLite format
          sim::ParticleList droppedParticleList = fparticleListAction->YieldDroppedList();
          droppedPartCol->reserve(droppedParticleList.size());

          for (auto const& partPair : droppedParticleList) {
            simb::MCParticle& p = *(partPair.second);
            if (ParticleListAction::isDropped(&p)) continue;
            if (p.StatusCode() != 1) continue;

            sim::MCParticleLite mini_mcp(p);
            mini_mcp.Origin(mct->Origin());

            droppedPartCol->push_back(std::move(mini_mcp));
          } // for(droppedParticleList)
        }

        // Has the user request a detailed dump of the output objects?
        if (fdumpParticleList) {
          mf::LogInfo("LArG4") << "Dump sim::ParticleList; size()=" << particleList.size() << "\n"
                               << particleList;
        }
      }

    } // end loop over interactions

    // get the electrons from the LArVoxelReadout sensitive detector
    // Get the sensitive-detector manager.
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    // Find the sensitive detector with the name "LArVoxelSD".
    auto theOpDetDet = dynamic_cast<OpDetSensitiveDetector*>(
      sdManager->FindSensitiveDetector("OpDetSensitiveDetector"));

    // Store the contents of the detected photon table
    //
    if (theOpDetDet) {

      if (!lgp->NoPhotonPropagation()) {

        for (int Reflected = 0; Reflected <= 1; Reflected++) {
          if (Reflected && !fStoreReflected) continue;

          if (!fUseLitePhotons) {
            MF_LOG_DEBUG("Optical") << "Storing OpDet Hit Collection in Event";
            std::vector<sim::SimPhotons>& ThePhotons =
              OpDetPhotonTable::Instance()->GetPhotons(Reflected);
            if (Reflected)
              PhotonColRefl->reserve(ThePhotons.size());
            else
              PhotonCol->reserve(ThePhotons.size());
            for (auto& it : ThePhotons) {
              if (Reflected)
                PhotonColRefl->push_back(std::move(it));
              else
                PhotonCol->push_back(std::move(it));
            }
          }
          else {
            MF_LOG_DEBUG("Optical") << "Storing OpDet Hit Collection in Event";

            std::map<int, std::map<int, int>> ThePhotons =
              OpDetPhotonTable::Instance()->GetLitePhotons(Reflected);

            if (size(ThePhotons) > 0) {
              LitePhotonCol->reserve(ThePhotons.size());
              for (auto const& [opChannel, detectedPhotons] : ThePhotons) {
                sim::SimPhotonsLite ph;
                ph.OpChannel = opChannel;
                ph.DetectedPhotons = detectedPhotons;
                if (Reflected)
                  LitePhotonColRefl->push_back(std::move(ph));
                else
                  LitePhotonCol->push_back(std::move(ph));
              }
            }
          }
          if (Reflected)
            *cOpDetBacktrackerRecordColRefl =
              OpDetPhotonTable::Instance()->YieldReflectedOpDetBacktrackerRecords();
          else
            *cOpDetBacktrackerRecordCol =
              OpDetPhotonTable::Instance()->YieldOpDetBacktrackerRecords();
        }
      } //end if no photon propagation

      if (lgp->FillSimEnergyDeposits()) {
        // we steal the only existing copy of the energy deposit map. Oink!
        auto edepMap = OpDetPhotonTable::Instance()->YieldSimEnergyDeposits();
        for (auto& [volumeName, edepCol] : edepMap) {
          // note: constant reference to a (smart) pointer to non-const data
          auto const& destColl =
            boost::contains(volumeName, "TPCActive") ? edepCol_TPCActive : edepCol_Other;
          append(*destColl, std::move(edepCol));
        } // for
      }
    } //end if theOpDetDet

    if (!lgp->NoElectronPropagation()) {

      // only put the sim::SimChannels into the event once, not once for every
      // MCTruth in the event

      std::set<LArVoxelReadout*> ReadoutList; // to be cleared later on

      for (unsigned int c = 0; c < geom->Ncryostats(); ++c) {

        // map to keep track of which channels we already have SimChannels for in scCol
        // remake this map on each cryostat as channels ought not to be shared between
        // cryostats, just between TPC's

        std::map<unsigned int, unsigned int> channelToscCol;

        unsigned int ntpcs = geom->Cryostat(c).NTPC();
        for (unsigned int t = 0; t < ntpcs; ++t) {
          std::string name("LArVoxelSD");
          std::ostringstream sstr;
          sstr << name << "_Cryostat" << c << "_TPC" << t;

          // try first to find the sensitive detector specific for this TPC;
          // do not bother writing on screen if there is none (yet)
          G4VSensitiveDetector* sd = sdManager->FindSensitiveDetector(sstr.str(), false);
          // if there is none, catch the general one (called just "LArVoxelSD")
          if (!sd) sd = sdManager->FindSensitiveDetector(name, false);
          // If this didn't work, then a sensitive detector with
          // the name "LArVoxelSD" does not exist.
          if (!sd) {
            throw cet::exception("LArG4")
              << "Sensitive detector for cryostat " << c << " TPC " << t << " not found (neither '"
              << sstr.str() << "' nor '" << name << "' exist)\n";
          }

          // Convert the G4VSensitiveDetector* to a LArVoxelReadout*.
          auto larVoxelReadout = dynamic_cast<LArVoxelReadout*>(sd);

          // If this didn't work, there is a "LArVoxelSD" detector, but
          // it's not a LArVoxelReadout object.
          if (!larVoxelReadout) {
            throw cet::exception("LArG4")
              << "Sensitive detector '" << sd->GetName() << "' is not a LArVoxelReadout object\n";
          }

          LArVoxelReadout::ChannelMap_t& channels = larVoxelReadout->GetSimChannelMap(c, t);
          if (!empty(channels)) {
            MF_LOG_DEBUG("LArG4") << "now put " << channels.size() << " SimChannels from C=" << c
                                  << " T=" << t << " into the event";
          }

          for (auto ch_pair : channels) {
            sim::SimChannel& sc = ch_pair.second;

            // push sc onto scCol but only if we haven't already put something in scCol for this channel.
            // if we have, then merge the ionization deposits.  Skip the check if we only have one TPC

            if (ntpcs > 1) {
              unsigned int ichan = sc.Channel();
              auto itertest = channelToscCol.find(ichan);
              if (itertest == channelToscCol.end()) {
                channelToscCol[ichan] = scCol->size();
                scCol->emplace_back(std::move(sc));
              }
              else {
                unsigned int idtest = itertest->second;
                auto const& tdcideMap = sc.TDCIDEMap();
                for (auto const& tdcide : tdcideMap) {
                  for (auto const& ide : tdcide.second) {
                    double xyz[3] = {ide.x, ide.y, ide.z};
                    scCol->at(idtest).AddIonizationElectrons(ide.trackID,
                                                             tdcide.first,
                                                             ide.numElectrons,
                                                             xyz,
                                                             ide.energy,
                                                             ide.origTrackID);
                  } // end loop to add ionization electrons to  scCol->at(idtest)
                }   // end loop over tdc to vector<sim::IDE> map
              }     // end if check to see if we've put SimChannels in for ichan yet or not
            }
            else {
              scCol->emplace_back(std::move(sc));
            } // end of check if we only have one TPC (skips check for multiple simchannels if we have just one TPC)
          }   // end loop over simchannels for this TPC

          // mark it for clearing
          ReadoutList.insert(const_cast<LArVoxelReadout*>(larVoxelReadout));

        } // end loop over tpcs
      }   // end loop over cryostats

      for (LArVoxelReadout* larVoxelReadout : ReadoutList) {
        larVoxelReadout->ClearSimChannels();
      }
    } //endif electron prop

    // only put the sim::AuxDetSimChannels into the event once, not once for every
    // MCTruth in the event

    adCol->reserve(geom->NAuxDets());
    for (unsigned int a = 0; a < geom->NAuxDets(); ++a) {

      // there should always be at least one senstive volume because
      // we make one for the full aux det if none are specified in the
      // gdml file - see AuxDetGeo.cxx
      for (size_t sv = 0; sv < geom->AuxDet(a).NSensitiveVolume(); ++sv) {

        // N.B. this name convention is used when creating the
        //      AuxDetReadout SD in AuxDetReadoutGeometry
        std::stringstream name;
        name << "AuxDetSD_AuxDet" << a << "_" << sv;
        G4VSensitiveDetector* sd = sdManager->FindSensitiveDetector(name.str().c_str());
        if (!sd) {
          throw cet::exception("LArG4")
            << "Sensitive detector '" << name.str() << "' does not exist\n";
        }

        // Convert the G4VSensitiveDetector* to a AuxDetReadout*.
        larg4::AuxDetReadout* auxDetReadout = dynamic_cast<larg4::AuxDetReadout*>(sd);

        MF_LOG_DEBUG("LArG4") << "now put the AuxDetSimTracks in the event";

        const sim::AuxDetSimChannel adsc = auxDetReadout->GetAuxDetSimChannel();
        adCol->push_back(adsc);
        auxDetReadout->clear();
      }
    } // Loop over AuxDets

    if (partCol) {
      mf::LogInfo("LArG4") << "Geant4 simulated " << nGeneratedParticles
                           << " MC particles, we keep " << partCol->size() << " .";
    }

    if (fdumpSimChannels) {
      mf::LogVerbatim("DumpSimChannels")
        << "Event " << evt.id() << ": " << scCol->size() << " channels with signal";
      unsigned int nChannels = 0;
      for (const sim::SimChannel& sc : *scCol) {
        mf::LogVerbatim out("DumpSimChannels");
        out << " #" << nChannels << ": ";
        // dump indenting with "    ", but not on the first line
        sc.Dump(out, "  ");
        ++nChannels;
      } // for
    }   // if dump SimChannels

    if (!lgp->NoElectronPropagation()) evt.put(std::move(scCol));

    evt.put(std::move(adCol));
    if (partCol) evt.put(std::move(partCol));
    if (droppedPartCol) {
      std::cout << "LArG4 dropped particles length = " << droppedPartCol->size() << std::endl;
      evt.put(std::move(droppedPartCol));
    }
    if (tpassn) evt.put(std::move(tpassn));
    if (!lgp->NoPhotonPropagation()) {
      if (!fUseLitePhotons) {
        evt.put(std::move(PhotonCol));
        if (fStoreReflected) evt.put(std::move(PhotonColRefl), "Reflected");
      }
      else {
        evt.put(std::move(LitePhotonCol));
        evt.put(std::move(cOpDetBacktrackerRecordCol));
        if (fStoreReflected) {
          evt.put(std::move(LitePhotonColRefl), "Reflected");
          evt.put(std::move(cOpDetBacktrackerRecordColRefl), "Reflected");
        }
      }
    }

    if (lgp->FillSimEnergyDeposits()) {
      evt.put(std::move(edepCol_TPCActive), "TPCActive");
      evt.put(std::move(edepCol_Other), "Other");
    }
    return;
  } // LArG4::produce()

} // namespace LArG4

DEFINE_ART_MODULE(larg4::LArG4)
