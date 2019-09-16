/**
 * @file   larsim/EventGenerator/LightSource_module.cc
 * @author Ben Jones, MIT 2010
 */

/**
 * @class evgen::LightSource
 * @brief Light source event generator which simulate an extended isotropic photon source
 *
 * The light source can be run in two modes, file mode or scan mode.  Each requires
 * the specification of a different set of parameters.
 *
 * File mode
 * ----------
 *
 * Light source position, intensity and shape are supplied on an event by event basis
 * in a text file.  See the example provided for the format. Pararmeters required:
 *
 *     int32   SourceMode = 0      - sets light source to file mode
 *     string  FileName            - file of per event light source specifications
 *     int32   PosDist             - how to distribute production points sampled in momentum, position
 *     int32   PDist                   and time ranges specified.  For all of these :
 *     int32   TDist                   0 = uniform and 1 = gauss
 *     bool    FillTree            - whether to write a tree of photon production points to fileservice
 *
 * Upon reaching the end of the file, the light source will loop back to the first point.
 * hence a one line text file will give a constant light source size, position and intensity.
 *
 * Scan mode
 * ----------
 *
 * Divide volume into cuboidal regions and produce an isotropic light source in each,
 * using one region per event.  User can specify either to use the full detector volume
 * or some custom specified volume.
 *
 * This mode is used when building a fast photon sim library, and performing volume
 * scan sensitivity studies.
 *
 *     int32   SourceMode = 1      - sets light source to scan mode
 *     int32   N                   - number of photons to shoot from each point
 *     double  P                   - peak photon momentum (or energy) in eV
 *     double  SigmaP              - momentum distribution width
 *     double  XSteps              - Number of regions to divide volume into in each direction
 *     double  YSteps
 *     double  ZSteps
 *     double  T0                  - Peak time of photon production
 *     double  SigmaT              - time distribution width
 *     int32   PosDist             - how to distribute production points sampled in momentum, position
 *     int32   PDist                 and time ranges specified.  For all of these :
 *     int32   TDist                   0 = uniform and 1 = gaussian
 *     bool    FillTree            - whether to write a tree of photon production points to fileservice
 *     bool    UseCustomRegion     - supply our own volme specification or use the full detector volume?
 *     vdouble[3]  RegionMin       - bounding corners of the custom volume specification
 *     vdouble[3]  RegionMax           (only used if UseCustomRegion=true)
 *
 */

// C++ includes.
#include <string>
#include <cmath>
#include <memory>
#include <fstream>

// ART includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/exception.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "larsim/Simulation/PhotonVoxels.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace evgen {

  /// A module for optical MC testing and library building
  class LightSource : public art::EDProducer {
  public:
    explicit LightSource(fhicl::ParameterSet const& pset);

    void produce(art::Event & evt);
    void beginRun(art::Run& run);

  private:

    void Sample(simb::MCTruth &truth);

    // for c2: fSeed is unused
    //int               fSeed;              //random number seed
    std::string       fVersion;           //version of the configuration

    // Flags to mark module modes
    static const int  kUNIF = 0;
    static const int  kGAUS = 1;
    static const int  kFILE = 0;
    static const int  kSCAN = 1;

    // File stream, filename and empty string for file processing
    std::ifstream      fInputFile;
    std::string        fFileName;
    char               fDummyString[256];

    // A ttree to keep track of where particles have been shot - ends up in histos.root
    TTree *            fPhotonsGenerated;
    TLorentzVector     fShotPos;
    TLorentzVector     fShotMom;
    Int_t              fEvID;

    // Parameters loaded from config - both modes
    int                fSourceMode;     // Mode to run in - scan or file
    bool               fFillTree;       // Do we want to create a TTree of shot particles?
    int                fPosDist;        //
    int                fTDist;          // Random distributions to use : 1= gauss, 0= uniform
    int                fPDist;          //

    //Scan mode specific parameters
    int fXSteps;                        //
    int fYSteps;                        //  Number of steps to take in each dimension
    int fZSteps;                        //

    sim::PhotonVoxelDef fThePhotonVoxelDef;  // The photon voxel definition object for scan mode

    int fVoxelCount;                    // Total number of voxels
    int fCurrentVoxel;                  // Counter to keep track of vox ID


    //  TPC Measurements
    TVector3 fTPCCenter;
    TVector3 fTPCDimensions;
    std::vector<double> fRegionMin;
    std::vector<double> fRegionMax;
    bool fUseCustomRegion;


    // Parameters used to shoot in distributions
    double              fX;              // central x position of source
    double              fY;              // central y position of source
    double              fZ;              // central z position of source
    bool                fPointSource;    // Point-like light source in fX, fY, fZ
    double              fT;              // central t position of source
    double              fSigmaX;         // x width
    double              fSigmaY;         // y width
    double              fSigmaZ;         // z width
    double              fSigmaT;         // t width
    double              fP;              // central momentm of photon
    double              fSigmaP;         // mom width;

    // Number of photons per event
    int                fN;              // number of photons per event

    int                fFirstVoxel;
    int                fLastVoxel;
    CLHEP::HepRandomEngine& fEngine;
  };
}

namespace evgen{

  //----------------------------------------------------------------
  LightSource::LightSource(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fSourceMode{pset.get<int >("SourceMode")}
    , fFillTree{pset.get<bool>("FillTree")}
    , fPosDist{pset.get<int >("PosDist")}
    , fTDist{pset.get<int >("TDist")}
    , fPDist{pset.get<int >("PDist")}
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    , fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, pset, "Seed"))
  {
    // load optional parameters in function
    produces< sumdata::RunData, art::InRun >();
    produces< std::vector<simb::MCTruth> >();

    if(fSourceMode==kFILE)
      {
	fFileName  = pset.get<std::string>("SteeringFile");
	fInputFile.open(fFileName.c_str());
	fInputFile.getline(fDummyString,256);

      }
    else if (fSourceMode==kSCAN)
      {
	fT      = pset.get<double>("T0");
	fSigmaT = pset.get<double>("SigmaT");
	fN      = pset.get<int   >("N");

	fFirstVoxel     = pset.get<int   >("FirstVoxel");
	fLastVoxel      = pset.get<int   >("LastVoxel");

	fP      = pset.get<double>("P");
	fSigmaP = pset.get<double>("SigmaP");

	fUseCustomRegion = pset.get<bool>("UseCustomRegion");
	fPointSource = pset.get<bool>("PointSource",false);

	if(fUseCustomRegion)
	  {
	    fRegionMin = pset.get< std::vector<double> >("RegionMin");
	    fRegionMax = pset.get< std::vector<double> >("RegionMax");
	    fXSteps = pset.get<int >("XSteps");
	    fYSteps = pset.get<int >("YSteps");
	    fZSteps = pset.get<int >("ZSteps");
	  }

        art::ServiceHandle<geo::Geometry const> geo;
	// get TPC dimensions removed. -TA


	fCurrentVoxel=0;

	// define voxelization based on parameters read from config.
	// There are two modes - either read the dimensions of the TPC from
	// the geometry, or use values specified by the user.
	if(!fUseCustomRegion)
	  {
            art::ServiceHandle<phot::PhotonVisibilityService const> vis;
	    fThePhotonVoxelDef = vis->GetVoxelDef();
	  }
	else
	  {
	    fThePhotonVoxelDef = sim::PhotonVoxelDef(fRegionMin[0],
						     fRegionMax[0],
						     fXSteps,
						     fRegionMin[1],
						     fRegionMax[1],
						     fYSteps,
						     fRegionMin[2],
						     fRegionMax[2],
						     fZSteps);
	  }


	// Set distribution widths to voxel size

	fSigmaX = fThePhotonVoxelDef.GetVoxelSize().X()/2.0;
	fSigmaY = fThePhotonVoxelDef.GetVoxelSize().Y()/2.0;
	fSigmaZ = fThePhotonVoxelDef.GetVoxelSize().Z()/2.0;

	// Get number of voxels we will step through

	fVoxelCount = fThePhotonVoxelDef.GetNVoxels();

	if(fLastVoxel<0) fLastVoxel = fVoxelCount;

	mf::LogVerbatim("LightSource") << "Light Source : Determining voxel params : "
				       << fVoxelCount << " "
				       << fSigmaX     << " "
				       << fSigmaY     << " "
				       <<fSigmaZ;

      }
    else{
      throw cet::exception("LightSource") << "EVGEN Light Source : Unrecognised light source mode\n";
    }

    if(fFillTree)
      {
        art::ServiceHandle<art::TFileService const> tfs;
	fPhotonsGenerated = tfs->make<TTree>("PhotonsGenerated","PhotonsGenerated");
	fPhotonsGenerated->Branch("X",&(fShotPos[0]),"X/D");
	fPhotonsGenerated->Branch("Y",&(fShotPos[1]),"Y/D");
	fPhotonsGenerated->Branch("Z",&(fShotPos[2]),"Z/D");
	fPhotonsGenerated->Branch("T",&(fShotPos[3]),"T/D");
	fPhotonsGenerated->Branch("PX",&(fShotMom[0]),"PX/D");
	fPhotonsGenerated->Branch("PY",&(fShotMom[1]),"PY/D");
	fPhotonsGenerated->Branch("PZ",&(fShotMom[2]),"PZ/D");
	fPhotonsGenerated->Branch("PT",&(fShotMom[3]),"PT/D");
	fPhotonsGenerated->Branch("EventID",&fEvID,"EventID/I");
      }
  }


  //____________________________________________________________________________
  void LightSource::beginRun(art::Run& run)
  {
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));

    fCurrentVoxel=fFirstVoxel;
  }

  //----------------------------------------------------------------
  void LightSource::produce(art::Event& evt)
  {
    if(fSourceMode==kFILE) {
    //  Each event, read coordinates of gun and number of photons to shoot from file
	// Loop file if required
	if(fInputFile.eof()){
	  mf::LogWarning("LightSource") << "EVGEN Light Source : Warning, reached end of file,"
					<< " looping back to beginning";
	  fInputFile.seekg(0,std::ios::beg);
	  fInputFile.clear();
	}

	if(!fInputFile.is_open() || fInputFile.fail() ){
	  throw cet::exception("LightSource") << "EVGEN Light Source : File error in "
					      << fFileName << "\n";
	}
	else{
	  // read in one line
	  fInputFile >> fX >> fY >> fZ >> fT
		     >> fSigmaX >> fSigmaY >> fSigmaZ >> fSigmaT
		     >> fP >> fSigmaP >> fN;
	  fInputFile.getline(fDummyString,256);
	  fThePhotonVoxelDef = sim::PhotonVoxelDef(fX - fSigmaX,
						   fX + fSigmaX,
						   1,
						   fY - fSigmaY,
						   fY + fSigmaY,
						   1,
						   fZ - fSigmaZ,
						   fZ + fSigmaZ,
						   1);

	  fCurrentVoxel=0;
	}
      }
    else if(fSourceMode==kSCAN) {
    //  Step through detector using a number of steps provided in the config file
    //  firing a constant number of photons from each point
	auto const VoxelCenter = fThePhotonVoxelDef.GetPhotonVoxel(fCurrentVoxel).GetCenter();
	fX = VoxelCenter.X();
	fY = VoxelCenter.Y();
	fZ = VoxelCenter.Z();
      }
    else{
      //  Neither file or scan mode, probably a config file error
      throw cet::exception("LightSource") <<"EVGEN : Light Source, unrecognised source mode\n";
    }

    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    Sample(truth);

    truthcol->push_back(truth);
    evt.put(std::move(truthcol));

    phot::PhotonVisibilityService* vis = nullptr;
    try {
      vis = art::ServiceHandle<phot::PhotonVisibilityService>().get();
    }
    catch (art::Exception const& e) {
      // if the service is not configured, then this is not a build job
      // (it is a simple generation job instead)
      if (e.categoryCode() != art::errors::ServiceNotFound) throw;
    }

    if(vis && vis->IsBuildJob())
      {
	mf::LogVerbatim("LightSource") << "Light source : Stowing voxel params ";
	vis->StoreLightProd(fCurrentVoxel,fN);
      }

    if(fCurrentVoxel!=fLastVoxel)
      {
	++fCurrentVoxel;
      }
    else
      {
	mf::LogVerbatim("LightSource") << "EVGEN Light Source fully scanned detector.  Starting over.";
	fCurrentVoxel=fFirstVoxel;
      }
  }


  void LightSource::Sample(simb::MCTruth& mct)
  {
    mf::LogVerbatim("LightSource") <<"Light source debug : Shooting at " << fX <<" " << fY<<" "<< fZ;

    CLHEP::RandFlat   flat(fEngine);
    CLHEP::RandGaussQ gauss(fEngine);

    for(int j=0; j!=fN; ++j){
      // Choose momentum (supplied in eV, convert to GeV)
      double p = fP;
      if (fPDist == kGAUS) {
	p = gauss.fire(fP, fSigmaP);
      }
      else {
	p = fP + fSigmaP*(2.0*flat.fire()-1.0);
      }
      p /= 1000000000.;

      // Choose position
      TVector3 x;
      if(fPointSource) {
	x[0] = fX;
	x[1] = fY;
	x[2] = fZ;
      }
      else {
        if (fPosDist == kGAUS) {
	  x[0] = gauss.fire(fX, fSigmaX);
	  x[1] = gauss.fire(fY, fSigmaY);
	  x[2] = gauss.fire(fZ, fSigmaZ);
        }
        else {
	  x[0] = fX + fSigmaX*(2.0*flat.fire()-1.0);
  	  x[1] = fY + fSigmaY*(2.0*flat.fire()-1.0);
  	  x[2] = fZ + fSigmaZ*(2.0*flat.fire()-1.0);
        }

      }

      // Choose time
      double t;
      if (fTDist == kGAUS) {
	t = gauss.fire(fT, fSigmaT);
      }
      else {
	t = fT + fSigmaT * (2.0 * flat.fire()-1.0);
      }


      //assume the position is relative to the center of the TPC
      //x += fTPCCenter;

      fShotPos = TLorentzVector(x[0], x[1], x[2], t);


      // Choose angles
      double costh = 2 * flat.fire() - 1;
      double sinth = pow(1-pow(costh,2),0.5);
      double phi   = 2 * M_PI * flat.fire();

      // Generate momentum 4-vector

      fShotMom = TLorentzVector( p*sinth*cos(phi),
				 p*sinth*sin(phi),
				 p*costh,
				 p                  );

      int trackid = -1*(j+1); // set track id to -i as these are all primary particles and have id <= 0
      std::string primary("primary");
      int PDG=0; //optical photons have PDG 0

      simb::MCParticle part(trackid, PDG, primary);
      part.AddTrajectoryPoint(fShotPos, fShotMom);

      if(fFillTree)
	fPhotonsGenerated->Fill();

      mct.Add(part);
    }

  }

}

DEFINE_ART_MODULE(evgen::LightSource)
