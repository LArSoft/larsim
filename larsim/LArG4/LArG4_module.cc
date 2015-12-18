////////////////////////////////////////////////////////////////////////
/// \file  LArG4.cxx
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \version $Id: LArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
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

#include "G4Base/G4Helper.h"
#include "G4Base/ConvertMCTruthToG4.h"

// C++ Includes
#include <sstream> // std::ostringstream
#include <vector> // std::ostringstream
#include <map> // std::ostringstream
#include <set> // std::ostringstream
#include <iostream>
// #include <cstring>
#include <sys/stat.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/Assns.h"
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

// art extensions
#include "artextensions/SeedService/SeedService.hh"

// LArSoft Includes
#include "LArG4/LArVoxelReadoutGeometry.h"
#include "LArG4/PhysicsList.h"
#include "LArG4/ParticleListAction.h"
#include "LArG4/G4BadIdeaAction.h"
#include "LArG4/IonizationAndScintillationAction.h"
#include "LArG4/OpDetSensitiveDetector.h"
#include "LArG4/OpDetReadoutGeometry.h"
#include "LArG4/LArStackingAction.h"
#include "LArG4/LArVoxelReadout.h"
#include "LArG4/MaterialPropertyLoader.h"
#include "LArG4/OpDetPhotonTable.h"
#include "LArG4/AuxDetReadoutGeometry.h"
#include "LArG4/AuxDetReadout.h"
#include "Simulation/LArG4Parameters.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationBase/MCTruth.h"
#include "Simulation/ParticleList.h"
#include "Simulation/SimPhotons.h"
#include "Simulation/SimChannel.h"
#include "Simulation/AuxDetSimChannel.h"
#include "Geometry/Geometry.h"
#include "G4Base/DetectorConstruction.h"
#include "G4Base/UserActionManager.h"

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


// Forward declarations
class G4RunManager;
class G4UImanager;
class G4VisExecutive;

///Geant4 interface
namespace larg4 {  
 
  // Forward declarations within namespace.
  class LArVoxelListAction;
  class ParticleListAction;
  
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

  private:

    g4b::G4Helper*             fG4Help;             ///< G4 interface object                                           
    larg4::LArVoxelListAction* flarVoxelListAction; ///< Geant4 user action to accumulate LAr voxel information.
    larg4::ParticleListAction* fparticleListAction; ///< Geant4 user action to particle information.                   

    std::string                fG4PhysListName;     ///< predefined physics list to use if not making a custom one
    std::string                fG4MacroPath;        ///< directory path for Geant4 macro file to be 
                                                    ///< executed before main MC processing.
    bool                       fdumpParticleList;   ///< Whether each event's sim::ParticleList will be displayed.
    bool                       fdumpSimChannels;    ///< Whether each event's sim::Channel will be displayed.
    bool                       fUseLitePhotons;
    int                        fSmartStacking;      ///< Whether to instantiate and use class to 
                                                    ///< dictate how tracks are put on stack.        
    std::vector<std::string>   fInputLabels;
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
    , fdumpParticleList      (pset.get< bool        >("DumpParticleList")                   )
    , fdumpSimChannels       (pset.get< bool        >("DumpSimChannels", false)             )
    , fSmartStacking         (pset.get< int         >("SmartStacking",0)                    )
  {
    LOG_DEBUG("LArG4") << "Debug: LArG4()";

    // setup the random number service for Geant4, the "G4Engine" label is a
    // special tag setting up a global engine for use by Geant4/CLHEP;
    // obtain the random seed from SeedService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<artext::SeedService>()
      ->createEngine(*this, "G4Engine", pset, "Seed");

    //get a list of generators to use, otherwise, we'll end up looking for anything that's
    //made an MCTruth object
    bool useInputLabels = pset.get_if_present< std::vector<std::string> >("InputLabels",fInputLabels);
    if(!useInputLabels) fInputLabels.resize(0);
    
    art::ServiceHandle<sim::LArG4Parameters> lgp;
    fUseLitePhotons = lgp->UseLitePhotons();

    if(!fUseLitePhotons) produces< std::vector<sim::SimPhotons>     >();
    else                 produces< std::vector<sim::SimPhotonsLite> >();

    produces< std::vector<simb::MCParticle> >();
    produces< std::vector<sim::SimChannel>  >();
    produces< std::vector<sim::AuxDetSimChannel> >();
    produces< art::Assns<simb::MCTruth, simb::MCParticle> >();

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

    fG4Help = new g4b::G4Helper(fG4MacroPath, fG4PhysListName);
    fG4Help->ConstructDetector(geom->GDMLFile());

    // Get the logical volume store and assign material properties
    larg4::MaterialPropertyLoader* MPL = new larg4::MaterialPropertyLoader();
    MPL->GetPropertiesFromServices();
    MPL->UpdateGeometry(G4LogicalVolumeStore::GetInstance());

    // Tell the detector about the parallel LAr voxel geometry.
    std::vector<G4VUserParallelWorld*> pworlds;

    // make a parallel world for each TPC in the detector
    pworlds.push_back( new LArVoxelReadoutGeometry("LArVoxelReadoutGeometry") );
    pworlds.push_back( new OpDetReadoutGeometry( geom->OpDetGeoName() ));
    pworlds.push_back( new AuxDetReadoutGeometry("AuxDetReadoutGeometry") );

    fG4Help->SetParallelWorlds(pworlds);

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

  //----------------------------------------------------------------------
  void LArG4::produce(art::Event& evt)
  {
    LOG_DEBUG("LArG4") << "produce()";

    // loop over the lists and put the particles and voxels into the event as collections
    std::unique_ptr< std::vector<simb::MCParticle> > partCol  (new std::vector<simb::MCParticle  >);
    std::unique_ptr< std::vector<sim::SimChannel>  > scCol    (new std::vector<sim::SimChannel>);
    std::unique_ptr< std::vector<sim::SimPhotons>  > PhotonCol(new std::vector<sim::SimPhotons>);
    std::unique_ptr< std::vector<sim::SimPhotonsLite>  > LitePhotonCol(new std::vector<sim::SimPhotonsLite>);
    std::unique_ptr< art::Assns<simb::MCTruth, simb::MCParticle> > tpassn(new art::Assns<simb::MCTruth, simb::MCParticle>);
    std::unique_ptr< std::vector< sim::AuxDetSimChannel > > adCol (new  std::vector<sim::AuxDetSimChannel> );

    // Fetch the lists of LAr voxels and particles.
    art::ServiceHandle<sim::LArG4Parameters> lgp;
    art::ServiceHandle<geo::Geometry> geom;

    // Clear the detected photon table
    OpDetPhotonTable::Instance()->ClearTable(geom->NOpDets());

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


    // Need to process Geant4 simulation for each interaction separately.
    for(size_t mcl = 0; mcl < mclists.size(); ++mcl){

      art::Handle< std::vector<simb::MCTruth> > mclistHandle = mclists[mcl];

      for(size_t m = 0; m < mclistHandle->size(); ++m){
        art::Ptr<simb::MCTruth> mct(mclistHandle, m);

        LOG_DEBUG("LArG4") << *(mct.get());

        // The following tells Geant4 to track the particles in this interaction.
        fG4Help->G4Run(mct);

        const sim::ParticleList& particleList = *( fparticleListAction->GetList() );
        
        partCol->reserve(partCol->size() + particleList.size()); // not very useful here...
        for(auto pitr = particleList.begin(); pitr != particleList.end(); ++pitr){
          // copy the particle so that it isnt const
          simb::MCParticle p(*(*pitr).second);
          partCol->push_back(p);
          util::CreateAssn(*this, evt, *(partCol.get()), mct, *(tpassn.get()));
        }

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
      }
    }
      
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
              scCol->push_back(std::move(sc));
            }
            else {
              unsigned int idtest = itertest->second;
              auto const tdcideMap = sc.TDCIDEMap();
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
            scCol->push_back(std::move(sc));
          } // end of check if we only have one TPC (skips check for multiple simchannels if we have just one TPC)
        } // end loop over simchannels for this TPC
        // mark it for clearing
        ReadoutList.insert(const_cast<LArVoxelReadout*>(larVoxelReadout));
      } // end loop over tpcs
    }// end loop over cryostats

    for (LArVoxelReadout* larVoxelReadout: ReadoutList)
      larVoxelReadout->ClearSimChannels();
    
    
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
					<< name
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
	
    if (fdumpSimChannels) {
      mf::LogInfo out("DumpSimChannels");
      out << "Wires: " << scCol->size() << " channels with signal" << std::endl;
      unsigned int nChannels = 0;
      for (const sim::SimChannel& sc: *scCol) {
        out << " #" << nChannels << ": ";
        // dump indenting with "    ", but not on the first line
        sc.Dump(out, "  ");
        ++nChannels;
      } // for
    } // if dump SimChannels
    evt.put(std::move(scCol));
    
    evt.put(std::move(adCol));
    evt.put(std::move(partCol));
    if(!fUseLitePhotons) evt.put(std::move(PhotonCol));
    else evt.put(std::move(LitePhotonCol));
    evt.put(std::move(tpassn));

    return;
  } // LArG4::produce()

} // namespace LArG4

namespace larg4 {

  DEFINE_ART_MODULE(LArG4)

} // namespace LArG4

#endif // LARG4_LARG4_H
