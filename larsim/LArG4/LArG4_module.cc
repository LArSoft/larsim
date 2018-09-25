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

#ifndef LARG4_LARG4_H
#define LARG4_LARG4_H 1

#include "nutools/G4Base/G4Helper.h"
#include "nutools/G4Base/ConvertMCTruthToG4.h"

// C++ Includes
#include <sstream> // std::ostringstream
#include <vector> // std::ostringstream
#include <map> // std::ostringstream
#include <set> // std::ostringstream
#include <iostream>
#include <sys/stat.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft Includes
#include "larsim/LArG4/LArVoxelReadoutGeometry.h"
#include "larsim/LArG4/PhysicsList.h"
#include "larsim/LArG4/ParticleListAction.h"
#include "larsim/LArG4/G4BadIdeaAction.h"
#include "larsim/LArG4/IonizationAndScintillationAction.h"
#include "larsim/LArG4/OpDetSensitiveDetector.h"
#include "larsim/LArG4/OpDetReadoutGeometry.h"
#include "larsim/LArG4/LArStackingAction.h"
#include "larsim/LArG4/LArVoxelReadout.h"
#include "larsim/LArG4/MaterialPropertyLoader.h"
#include "larsim/LArG4/OpDetPhotonTable.h"
#include "larsim/LArG4/AuxDetReadoutGeometry.h"
#include "larsim/LArG4/AuxDetReadout.h"
#include "larsim/LArG4/ParticleFilters.h" // larg4::PositionInVolumeFilter
#include "lardataalg/MCDumpers/MCDumpers.h" // sim::dump namespace
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataobj/Simulation/sim.h"
#include "larcore/Geometry/Geometry.h"
#include "nutools/G4Base/DetectorConstruction.h"
#include "nutools/G4Base/UserActionManager.h"

// G4 Includes
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4VUserDetectorConstruction.hh"
#include "Geant4/G4VUserPrimaryGeneratorAction.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4UserRunAction.hh"
#include "Geant4/G4UserEventAction.hh"
#include "Geant4/G4UserTrackingAction.hh"
#include "Geant4/G4UserSteppingAction.hh"
#include "Geant4/G4UserStackingAction.hh"
#include "Geant4/G4VisExecutive.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4LogicalVolumeStore.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/globals.hh"

// ROOT Includes
#include "TGeoManager.h"

//For energy depositions
#include "lardataobj/Simulation/SimEnergyDeposit.h"


// Forward declarations
class G4RunManager;
class G4UImanager;
class G4VisExecutive;

#if defined __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-private-field"
#endif

///Geant4 interface
namespace larg4 {  
 
  // Forward declarations within namespace.
  class LArVoxelListAction;
  class ParticleListAction;
  
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
   * The interface with Geant4 is via a helper class provided by _nutools_.
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
   *   The interface with Geant4 is via a helper class provided by _nutools_.
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
   *     whether to print all MCParticles tracked
   * - *DumpSimChannels* (bool, default: `false`):
   *     whether to print all depositions on each SimChannel
   * - *SmartStacking* (int, default: `0`):
   *     whether to use class to dictate how tracks are put on stack (nonzero is on)
   * - *KeepParticlesInVolumes* (list of strings, default: _empty_):
   *     list of volumes in which to keep `simb::MCParticle` objects (empty keeps all)
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
   */
  class LArG4 : public art::EDProducer{
  public:
 
    /// Standard constructor and destructor for an FMWK module.
    explicit LArG4(fhicl::ParameterSet const& pset);
    virtual ~LArG4();

    /// The main routine of this module: Fetch the primary particles
    /// from the event, simulate their evolution in the detctor, and
    /// produce the detector response.
    void produce (art::Event& evt); 
    void beginJob();
    void beginRun(art::Run& run);

  private:
    g4b::G4Helper*             fG4Help;             ///< G4 interface object                                           
    larg4::LArVoxelListAction* flarVoxelListAction; ///< Geant4 user action to accumulate LAr voxel information.
    larg4::ParticleListAction* fparticleListAction; ///< Geant4 user action to particle information.                   

    std::string                fG4PhysListName;     ///< predefined physics list to use if not making a custom one
    std::string                fG4MacroPath;        ///< directory path for Geant4 macro file to be 
                                                    ///< executed before main MC processing.
    bool                       fCheckOverlaps;      ///< Whether to use the G4 overlap checker
    bool                       fdumpParticleList;   ///< Whether each event's sim::ParticleList will be displayed.
    bool                       fdumpSimChannels;    ///< Whether each event's sim::Channel will be displayed.
    bool                       fUseLitePhotons;
    int                        fSmartStacking;      ///< Whether to instantiate and use class to 
    double                     fOffPlaneMargin = 0.; ///< Off-plane charge recovery margin
                                                    ///< dictate how tracks are put on stack.        
    std::vector<std::string>   fInputLabels;
    std::vector<std::string>   fKeepParticlesInVolumes; ///<Only write particles that have trajectories through these volumes
    
    bool fSparsifyTrajectories; ///< Sparsify MCParticle Trajectories

    /// Configures and returns a particle filter
    std::unique_ptr<PositionInVolumeFilter> CreateParticleVolumeFilter
      (std::set<std::string> const& vol_names) const;
    
  };

} // namespace LArG4

namespace larg4 {

  //----------------------------------------------------------------------
  // Constructor
  LArG4::LArG4(fhicl::ParameterSet const& pset)
    : fG4Help                (0)
    , flarVoxelListAction    (0)
    , fparticleListAction    (0)
    , fG4PhysListName        (pset.get< std::string >("G4PhysListName","larg4::PhysicsList"))
    , fCheckOverlaps         (pset.get< bool        >("CheckOverlaps",false)                )
    , fdumpParticleList      (pset.get< bool        >("DumpParticleList",false)             )
    , fdumpSimChannels       (pset.get< bool        >("DumpSimChannels", false)             )
    , fSmartStacking         (pset.get< int         >("SmartStacking",0)                    )
    , fOffPlaneMargin        (pset.get< double      >("ChargeRecoveryMargin",0.0)           )
    , fKeepParticlesInVolumes        (pset.get< std::vector< std::string > >("KeepParticlesInVolumes",{}))
    , fSparsifyTrajectories  (pset.get< bool        >("SparsifyTrajectories",false)         )

  {
    LOG_DEBUG("LArG4") << "Debug: LArG4()";
    art::ServiceHandle<art::RandomNumberGenerator> rng;

    if (pset.has_key("Seed")) {
      throw art::Exception(art::errors::Configuration)
        << "The configuration of LArG4 module has the discontinued 'Seed' parameter.\n"
        "Seeds are now controlled by two parameters: 'GEANTSeed' and 'PropagationSeed'.";
    }
    // setup the random number service for Geant4, the "G4Engine" label is a
    // special tag setting up a global engine for use by Geant4/CLHEP;
    // obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed" or "GEANTSeed"
    art::ServiceHandle<rndm::NuRandomService>()
      ->createEngine(*this, "G4Engine", "GEANT", pset, "GEANTSeed");
    // same thing for the propagation engine:
    art::ServiceHandle<rndm::NuRandomService>()
      ->createEngine(*this, "HepJamesRandom", "propagation", pset, "PropagationSeed");

    //get a list of generators to use, otherwise, we'll end up looking for anything that's
    //made an MCTruth object
    bool useInputLabels = pset.get_if_present< std::vector<std::string> >("InputLabels",fInputLabels);
    if(!useInputLabels) fInputLabels.resize(0);
    
    art::ServiceHandle<sim::LArG4Parameters> lgp;

    fUseLitePhotons = lgp->UseLitePhotons();
    if(!lgp->NoPhotonPropagation()){
      if(!fUseLitePhotons) produces< std::vector<sim::SimPhotons>     >();
      else{
	produces< std::vector<sim::SimPhotonsLite> >();
	produces< std::vector<sim::OpDetBacktrackerRecord>   >();
      }
    }

    if(lgp->FillSimEnergyDeposits()){
      produces < std::vector<sim::SimEnergyDeposit> >("TPCActive");
      produces < std::vector<sim::SimEnergyDeposit> >("Other");
    }

    produces< std::vector<simb::MCParticle> >();
    if(!lgp->NoElectronPropagation()) produces< std::vector<sim::SimChannel>  >();
    produces< std::vector<sim::AuxDetSimChannel> >();
    produces< art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo> >();

    // constructor decides if initialized value is a path or an environment variable
    cet::search_path sp("FW_SEARCH_PATH");

    sp.find_file(pset.get< std::string >("GeantCommandFile"), fG4MacroPath);
    struct stat sb;
    if (fG4MacroPath.empty() || stat(fG4MacroPath.c_str(), &sb)!=0)
      // failed to resolve the file name
      throw cet::exception("NoG4Macro") << "G4 macro file "
                                        << fG4MacroPath
                                        << " not found!\n";

  }

  //----------------------------------------------------------------------
  // Destructor
  LArG4::~LArG4()
  {
    if(fG4Help) delete fG4Help;
  }

  //----------------------------------------------------------------------
  void LArG4::beginJob()
  {
    art::ServiceHandle<geo::Geometry> geom;
    auto* rng = &*(art::ServiceHandle<art::RandomNumberGenerator>());

    fG4Help = new g4b::G4Helper(fG4MacroPath, fG4PhysListName);
    if(fCheckOverlaps) fG4Help->SetOverlapCheck(true);
    fG4Help->ConstructDetector(geom->GDMLFile());

    // Get the logical volume store and assign material properties
    larg4::MaterialPropertyLoader* MPL = new larg4::MaterialPropertyLoader();
    MPL->GetPropertiesFromServices();
    MPL->UpdateGeometry(G4LogicalVolumeStore::GetInstance());

    // Tell the detector about the parallel LAr voxel geometry.
    std::vector<G4VUserParallelWorld*> pworlds;
    // Intialize G4 physics and primary generator action
    fG4Help->InitPhysics();

    // create the ionization and scintillation calculator;
    // this is a singleton (!) so it does not make sense
    // to create it in LArVoxelReadoutGeometry
    IonizationAndScintillation::CreateInstance(rng->getEngine("propagation"));

    // make a parallel world for each TPC in the detector
    LArVoxelReadoutGeometry::Setup_t readoutGeomSetupData;
    readoutGeomSetupData.readoutSetup.offPlaneMargin = fOffPlaneMargin;
    readoutGeomSetupData.readoutSetup.propGen
      = &(rng->getEngine("propagation"));
    pworlds.push_back(new LArVoxelReadoutGeometry
      ("LArVoxelReadoutGeometry", readoutGeomSetupData)
      );
    pworlds.push_back( new OpDetReadoutGeometry( geom->OpDetGeoName() ));
    pworlds.push_back( new AuxDetReadoutGeometry("AuxDetReadoutGeometry") );

    fG4Help->SetParallelWorlds(pworlds);

    // moved up
    // Intialize G4 physics and primary generator action
       fG4Help->InitPhysics();

    // Use the UserActionManager to handle all the Geant4 user hooks.
    g4b::UserActionManager* uaManager = g4b::UserActionManager::Instance();

    // User-action class for accumulating LAr voxels.
    art::ServiceHandle<sim::LArG4Parameters> lgp;

    // UserAction for getting past a bug in v4.9.4.p02 of Geant4.
    // This action will not be used once the bug has been fixed
    // The techniques used in this UserAction are not to be repeated
    // as in general they are a very bad idea, ie they take a const
    // pointer and jump through hoops to change it
    // 08-Apr-2014 WGS: It appears that with the shift to Geant 4.9.6 or
    // above, there's no longer any need for the "Bad Idea Action" fix.
    //    larg4::G4BadIdeaAction *bia = new larg4::G4BadIdeaAction(fSmartStacking);
    //    uaManager->AddAndAdoptAction(bia);

    // remove IonizationAndScintillationAction for now as we are ensuring
    // the Reset for each G4Step within the G4SensitiveVolumes
    //larg4::IonizationAndScintillationAction *iasa = new larg4::IonizationAndScintillationAction();
    //uaManager->AddAndAdoptAction(iasa);

    // User-action class for accumulating particles and trajectories
    // produced in the detector.
    fparticleListAction = new larg4::ParticleListAction(lgp->ParticleKineticEnergyCut(),
                                                        lgp->StoreTrajectories(),
                                                        lgp->KeepEMShowerDaughters());
    uaManager->AddAndAdoptAction(fparticleListAction);

    // UserActionManager is now configured so continue G4 initialization
    fG4Help->SetUserAction();

    // With an enormous detector with lots of rock ala LAr34 (nee LAr20)
    // we need to be smarter about stacking.
    if (fSmartStacking>0){
      G4UserStackingAction* stacking_action = new LArStackingAction(fSmartStacking);
      fG4Help->GetRunManager()->SetUserAction(stacking_action);
    }
    

  
  }

  void LArG4::beginRun(art::Run& run){
    // prepare the filter object (null if no filtering)
    
    std::set<std::string> volnameset(fKeepParticlesInVolumes.begin(), fKeepParticlesInVolumes.end());
    fparticleListAction->ParticleFilter(CreateParticleVolumeFilter(volnameset));
    
  }
  
  std::unique_ptr<PositionInVolumeFilter> LArG4::CreateParticleVolumeFilter
    (std::set<std::string> const& vol_names) const
  {
    
    // if we don't have favourite volumes, don't even bother creating a filter
    if (vol_names.empty()) return {};
    
    auto const& geom = *art::ServiceHandle<geo::Geometry>();
    
    std::vector<std::vector<TGeoNode const*>> node_paths
      = geom.FindAllVolumePaths(vol_names);
    
    // collection of interesting volumes
    PositionInVolumeFilter::AllVolumeInfo_t GeoVolumePairs;
    GeoVolumePairs.reserve(node_paths.size()); // because we are obsessed
  
    //for each interesting volume, follow the node path and collect
    //total rotations and translations
    for (size_t iVolume = 0; iVolume < node_paths.size(); ++iVolume) {
      std::vector<TGeoNode const*> path = node_paths[iVolume];
      
      TGeoTranslation* pTransl = new TGeoTranslation(0.,0.,0.);
      TGeoRotation* pRot = new TGeoRotation();
      for (TGeoNode const* node: path) {
        TGeoTranslation thistranslate(*node->GetMatrix());
        TGeoRotation thisrotate(*node->GetMatrix());
        pTransl->Add(&thistranslate);
        *pRot=*pRot * thisrotate;
      }
      
      //for some reason, pRot and pTransl don't have tr and rot bits set correctly
      //make new translations and rotations so bits are set correctly
      TGeoTranslation* pTransl2 = new TGeoTranslation(pTransl->GetTranslation()[0],
							pTransl->GetTranslation()[1],
						      pTransl->GetTranslation()[2]);
      double phi=0.,theta=0.,psi=0.;
      pRot->GetAngles(phi,theta,psi);
      TGeoRotation* pRot2 = new TGeoRotation();
      pRot2->SetAngles(phi,theta,psi);
      
      TGeoCombiTrans* pTransf = new TGeoCombiTrans(*pTransl2,*pRot2);

      GeoVolumePairs.emplace_back(node_paths[iVolume].back()->GetVolume(), pTransf);

    }
    
    return std::make_unique<PositionInVolumeFilter>(std::move(GeoVolumePairs));
    
  } // CreateParticleVolumeFilter()
    

  void LArG4::produce(art::Event& evt)
  {
    LOG_DEBUG("LArG4") << "produce()";

    // loop over the lists and put the particles and voxels into the event as collections
    std::unique_ptr< std::vector<sim::SimChannel>  >               scCol                      (new std::vector<sim::SimChannel>);
    std::unique_ptr< std::vector< sim::AuxDetSimChannel > >        adCol                      (new  std::vector<sim::AuxDetSimChannel> );
    auto tpassn = std::make_unique<art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo>>();
    std::unique_ptr< std::vector<simb::MCParticle> >               partCol                    (new std::vector<simb::MCParticle  >);
    std::unique_ptr< std::vector<sim::SimPhotons>  >               PhotonCol                  (new std::vector<sim::SimPhotons>);
    std::unique_ptr< std::vector<sim::SimPhotonsLite>  >           LitePhotonCol              (new std::vector<sim::SimPhotonsLite>);
    std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > >  cOpDetBacktrackerRecordCol (new std::vector<sim::OpDetBacktrackerRecord>);
    
    art::PtrMaker<simb::MCParticle> makeMCPartPtr(evt, *this);

    //for energy deposits
    std::unique_ptr< std::vector<sim::SimEnergyDeposit> > edepCol_TPCActive (new std::vector<sim::SimEnergyDeposit>);
    std::unique_ptr< std::vector<sim::SimEnergyDeposit> > edepCol_Other     (new std::vector<sim::SimEnergyDeposit>);

    // Fetch the lists of LAr voxels and particles.
    art::ServiceHandle<sim::LArG4Parameters> lgp;
    art::ServiceHandle<geo::Geometry> geom;

    // Clear the detected photon table
    OpDetPhotonTable::Instance()->ClearTable(geom->NOpDets());
    if(lgp->FillSimEnergyDeposits())
      OpDetPhotonTable::Instance()->ClearEnergyDeposits();

    // reset the track ID offset as we have a new collection of interactions
    fparticleListAction->ResetTrackIDOffset();

    //look to see if there is any MCTruth information for this
    //event
    std::vector< art::Handle< std::vector<simb::MCTruth> > > mclists;
    if(fInputLabels.size()==0)
      evt.getManyByType(mclists);
    else{
      mclists.resize(fInputLabels.size());
      for(size_t i=0; i<fInputLabels.size(); i++)
        evt.getByLabel(fInputLabels[i],mclists[i]);
    }
    
    unsigned int nGeneratedParticles = 0;
    
    // Need to process Geant4 simulation for each interaction separately.
    for(size_t mcl = 0; mcl < mclists.size(); ++mcl){

      art::Handle< std::vector<simb::MCTruth> > mclistHandle = mclists[mcl];

      for(size_t m = 0; m < mclistHandle->size(); ++m){
        art::Ptr<simb::MCTruth> mct(mclistHandle, m);

        LOG_DEBUG("LArG4") << *(mct.get());

        // The following tells Geant4 to track the particles in this interaction.
        fG4Help->G4Run(mct);

        // receive the particle list
        sim::ParticleList particleList = fparticleListAction->YieldList();
        
        for(auto const& partPair: particleList) {
          simb::MCParticle& p = *(partPair.second);
          ++nGeneratedParticles;
          
          // if the particle has been marked as dropped, we don't save it
          // (as of LArSoft ~v5.6 this does not ever happen because
          // ParticleListAction has already taken care of deleting them)
          if (ParticleListAction::isDropped(&p)) continue;
          
          sim::GeneratedParticleInfo const truthInfo{
            fparticleListAction->GetPrimaryTruthIndex(p.TrackId())
            };
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

	  if(fSparsifyTrajectories) p.SparsifyTrajectory();
          
          partCol->push_back(std::move(p));
          
          tpassn->addSingle(mct, makeMCPartPtr(partCol->size() - 1), truthInfo);
          
        } // for(particleList)


        // Has the user request a detailed dump of the output objects?
        if (fdumpParticleList){
          mf::LogInfo("LArG4") << "Dump sim::ParticleList; size()="
                               << particleList.size() << "\n"
                               << particleList;
        }

      }

    }// end loop over interactions
   
    // get the electrons from the LArVoxelReadout sensitive detector
    // Get the sensitive-detector manager.
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    
    // Find the sensitive detector with the name "LArVoxelSD".
    OpDetSensitiveDetector *theOpDetDet = dynamic_cast<OpDetSensitiveDetector*>(sdManager->FindSensitiveDetector("OpDetSensitiveDetector"));
 
    // Store the contents of the detected photon table
    //
    if(theOpDetDet){
      
      if(!lgp->NoPhotonPropagation()){

	if(!fUseLitePhotons){      
	  LOG_DEBUG("Optical") << "Storing OpDet Hit Collection in Event";
	  std::vector<sim::SimPhotons>& ThePhotons = OpDetPhotonTable::Instance()->GetPhotons();
	  PhotonCol->reserve(ThePhotons.size());
	  for(auto& it : ThePhotons)
	    PhotonCol->push_back(std::move(it));
	}
	else{
	  LOG_DEBUG("Optical") << "Storing OpDet Hit Collection in Event";
	  
	  std::map<int, std::map<int, int> > ThePhotons = OpDetPhotonTable::Instance()->GetLitePhotons();
	  
	  if(ThePhotons.size() > 0){
	    LitePhotonCol->reserve(ThePhotons.size());
	    for(auto const& it : ThePhotons){
	      sim::SimPhotonsLite ph;
	      ph.OpChannel = it.first;
	      ph.DetectedPhotons = it.second;
	      LitePhotonCol->push_back(ph);
	    }
	  }
	  *cOpDetBacktrackerRecordCol = OpDetPhotonTable::Instance()->YieldOpDetBacktrackerRecords();
	}
      }//end if no photon propagation

      if(lgp->FillSimEnergyDeposits())
	{
	  auto const& edepMap = OpDetPhotonTable::Instance()->GetSimEnergyDeposits();
	  for(auto const& edepCol : edepMap){
	    if(boost::contains(edepCol.first,"TPCActive"))
	      edepCol_TPCActive->insert(edepCol_TPCActive->end(),
					edepCol.second.begin(),edepCol.second.end());
	    else
	      edepCol_Other->insert(edepCol_Other->end(),
				    edepCol.second.begin(),edepCol.second.end());
	  }
	}
    }//end if theOpDetDet
    

    if(!lgp->NoElectronPropagation())
      {
    
	// only put the sim::SimChannels into the event once, not once for every
	// MCTruth in the event
	
	std::set<LArVoxelReadout*> ReadoutList; // to be cleared later on
	
	for(unsigned int c = 0; c < geom->Ncryostats(); ++c){
	  
	  // map to keep track of which channels we already have SimChannels for in scCol
	  // remake this map on each cryostat as channels ought not to be shared between 
	  // cryostats, just between TPC's
	  
	  std::map<unsigned int, unsigned int>  channelToscCol;
	  
	  unsigned int ntpcs =  geom->Cryostat(c).NTPC();
	  for(unsigned int t = 0; t < ntpcs; ++t){
	    std::string name("LArVoxelSD");
	    std::ostringstream sstr;
	    sstr << name << "_Cryostat" << c << "_TPC" << t;
	    
	    // try first to find the sensitive detector specific for this TPC;
	    // do not bother writing on screen if there is none (yet)
	    G4VSensitiveDetector* sd
	      = sdManager->FindSensitiveDetector(sstr.str(), false);
	    // if there is none, catch the general one (called just "LArVoxelSD")
	    if (!sd) sd = sdManager->FindSensitiveDetector(name, false);
	    // If this didn't work, then a sensitive detector with
	    // the name "LArVoxelSD" does not exist.
	    if ( !sd ){
	      throw cet::exception("LArG4") << "Sensitive detector for cryostat "
					    << c << " TPC " << t << " not found (neither '"
					    << sstr.str() << "' nor '" << name  << "' exist)\n";
	    }
	    
	    // Convert the G4VSensitiveDetector* to a LArVoxelReadout*.
	    LArVoxelReadout* larVoxelReadout = dynamic_cast<LArVoxelReadout*>(sd);
	    
	    // If this didn't work, there is a "LArVoxelSD" detector, but
	    // it's not a LArVoxelReadout object.
	    if ( !larVoxelReadout ){
	      throw cet::exception("LArG4") << "Sensitive detector '"
					    << sd->GetName()
					    << "' is not a LArVoxelReadout object\n";
	    }
	    
	    LArVoxelReadout::ChannelMap_t& channels = larVoxelReadout->GetSimChannelMap(c, t);
	    if (!channels.empty()) {
	      LOG_DEBUG("LArG4") << "now put " << channels.size() << " SimChannels"
		" from C=" << c << " T=" << t << " into the event";
	    }
	    
	    for(auto ch_pair: channels){
	      sim::SimChannel& sc = ch_pair.second;
	      
	      // push sc onto scCol but only if we haven't already put something in scCol for this channel.
	      // if we have, then merge the ionization deposits.  Skip the check if we only have one TPC
	      
	      if (ntpcs > 1) {
		unsigned int ichan = sc.Channel();
		std::map<unsigned int, unsigned int>::iterator itertest = channelToscCol.find(ichan);
		if (itertest == channelToscCol.end()) {
		  channelToscCol[ichan] = scCol->size();
		  scCol->emplace_back(std::move(sc));
		}
		else {
              unsigned int idtest = itertest->second;
              auto const& tdcideMap = sc.TDCIDEMap();
              for(auto const& tdcide : tdcideMap){
		for(auto const& ide : tdcide.second){
		  double xyz[3] = {ide.x, ide.y, ide.z};
		  scCol->at(idtest).AddIonizationElectrons(ide.trackID,
							   tdcide.first,
							   ide.numElectrons,
							   xyz,
							   ide.energy);
		} // end loop to add ionization electrons to  scCol->at(idtest)
	      }// end loop over tdc to vector<sim::IDE> map
		} // end if check to see if we've put SimChannels in for ichan yet or not
	      }
	      else {
		scCol->emplace_back(std::move(sc));
          } // end of check if we only have one TPC (skips check for multiple simchannels if we have just one TPC)
	    } // end loop over simchannels for this TPC
	    
	    
	    // mark it for clearing
	    ReadoutList.insert(const_cast<LArVoxelReadout*>(larVoxelReadout));
	    
	  } // end loop over tpcs
	}// end loop over cryostats
	
	for (LArVoxelReadout* larVoxelReadout: ReadoutList){
	  larVoxelReadout->ClearSimChannels();
	}
      }//endif electron prop
    
    // only put the sim::AuxDetSimChannels into the event once, not once for every
    // MCTruth in the event
    
    adCol->reserve(geom->NAuxDets());
    for(unsigned int a = 0; a < geom->NAuxDets(); ++a){

      // there should always be at least one senstive volume because 
      // we make one for the full aux det if none are specified in the 
      // gdml file - see AuxDetGeo.cxx
      for(size_t sv = 0; sv < geom->AuxDet(a).NSensitiveVolume(); ++sv){

          // N.B. this name convention is used when creating the
          //      AuxDetReadout SD in AuxDetReadoutGeometry
        std::stringstream name;
        name << "AuxDetSD_AuxDet" << a << "_" << sv;
        G4VSensitiveDetector* sd = sdManager->FindSensitiveDetector(name.str().c_str());
        if ( !sd ){
          throw cet::exception("LArG4") << "Sensitive detector '"
          << name.str()
          << "' does not exist\n";
        }
        
          // Convert the G4VSensitiveDetector* to a AuxDetReadout*.
        larg4::AuxDetReadout *auxDetReadout = dynamic_cast<larg4::AuxDetReadout*>(sd);
        
        LOG_DEBUG("LArG4") << "now put the AuxDetSimTracks in the event";
        
        const sim::AuxDetSimChannel adsc = auxDetReadout->GetAuxDetSimChannel();
        adCol->push_back(adsc);
        auxDetReadout->clear();
      }
      
    } // Loop over AuxDets
	
    mf::LogInfo("LArG4")
      << "Geant4 simulated " << nGeneratedParticles << " MC particles, we keep "
      << partCol->size() << " .";
    
    if (fdumpSimChannels) {
      mf::LogVerbatim("DumpSimChannels")
        << "Event " << evt.id()
        << ": " << scCol->size() << " channels with signal";
      unsigned int nChannels = 0;
      for (const sim::SimChannel& sc: *scCol) {
         mf::LogVerbatim out("DumpSimChannels");
         out << " #" << nChannels << ": ";
        // dump indenting with "    ", but not on the first line
        sc.Dump(out, "  ");
        ++nChannels;
      } // for
    } // if dump SimChannels

    if(!lgp->NoElectronPropagation()) evt.put(std::move(scCol));
    
    evt.put(std::move(adCol));
    evt.put(std::move(partCol));
    if(!lgp->NoPhotonPropagation()){
      if(!fUseLitePhotons) evt.put(std::move(PhotonCol));
      else{
	evt.put(std::move(LitePhotonCol));
	evt.put(std::move(cOpDetBacktrackerRecordCol));
      }
    }
    evt.put(std::move(tpassn));

    if(lgp->FillSimEnergyDeposits()){
      evt.put(std::move(edepCol_TPCActive),"TPCActive");
      evt.put(std::move(edepCol_Other),"Other");
    }
    return;
  } // LArG4::produce()

} // namespace LArG4

namespace larg4 {

  DEFINE_ART_MODULE(LArG4)

} // namespace LArG4
#if defined __clang__
  #pragma clang diagnostic pop
#endif

#endif // LARG4_LARG4_H
