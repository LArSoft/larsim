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
 *     bool    UseCustomRegion     - supply our own volume specification or use the full detector volume?
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
#include "cetlib/pow.h" // cet::square()
#include "cetlib_except/exception.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// nusimdata includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// lar includes
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "larsim/Simulation/PhotonVoxels.h"

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

    simb::MCTruth Sample();

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
    geo::Vector_t fTPCCenter;
    std::vector<double> fRegionMin;
    std::vector<double> fRegionMax;
    bool fUseCustomRegion;


    // Parameters used to shoot in distributions
    geo::Point_t        fCenter;         ///< Central position of source [cm]
    bool                fPointSource;    // Point-like light source in fCenter
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
          double x, y, z;
          fInputFile >> x >> y >> z >> fT
                     >> fSigmaX >> fSigmaY >> fSigmaZ >> fSigmaT
                     >> fP >> fSigmaP >> fN;
          fInputFile.getline(fDummyString,256);
          fCenter = { x, y, z };
          fThePhotonVoxelDef = sim::PhotonVoxelDef(fCenter.X() - fSigmaX,
                                                   fCenter.X() + fSigmaX,
                                                   1,
                                                   fCenter.Y() - fSigmaY,
                                                   fCenter.Y() + fSigmaY,
                                                   1,
                                                   fCenter.Z() - fSigmaZ,
                                                   fCenter.Z() + fSigmaZ,
                                                   1);

          fCurrentVoxel=0;
        }
      }
    else if(fSourceMode==kSCAN) {
    //  Step through detector using a number of steps provided in the config file
    //  firing a constant number of photons from each point
        fCenter = fThePhotonVoxelDef.GetPhotonVoxel(fCurrentVoxel).GetCenter();
      }
    else{
      //  Neither file or scan mode, probably a config file error
      throw cet::exception("LightSource") <<"EVGEN : Light Source, unrecognised source mode\n";
    }

    auto truthcol = std::make_unique<std::vector<simb::MCTruth>>();

    truthcol->push_back(Sample());
    
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


  simb::MCTruth LightSource::Sample()
  {
    mf::LogVerbatim("LightSource") <<"Light source debug : Shooting at " << fCenter;

    CLHEP::RandFlat   flat(fEngine, -1.0, 1.0);
    CLHEP::RandGaussQ gauss(fEngine);

    simb::MCTruth mct;
    mct.SetOrigin(simb::kSingleParticle);
    for(int j=0; j!=fN; ++j){
      // Choose momentum (supplied in eV, convert to GeV)
      double const p = 1e-9 * (
        (fPDist == kGAUS)
          ? gauss.fire(fP, fSigmaP)
          : fP + fSigmaP * flat.fire()
        );
      
      // Choose position
      geo::Point_t x;
      
      if(fPointSource) {
        x = fCenter;
      }
      else {
        if (fPosDist == kGAUS) {
          x = {
            gauss.fire(fCenter.X(), fSigmaX),
            gauss.fire(fCenter.Y(), fSigmaY),
            gauss.fire(fCenter.Z(), fSigmaZ)
          };
        }
        else {
          x = {
            fCenter.X() + fSigmaX * flat.fire(),
            fCenter.Y() + fSigmaY * flat.fire(),
            fCenter.Z() + fSigmaZ * flat.fire()
          };
        }

      }

      // Choose time
      double t;
      if (fTDist == kGAUS) {
        t = gauss.fire(fT, fSigmaT);
      }
      else {
        t = fT + fSigmaT * flat.fire();
      }


      //assume the position is relative to the center of the TPC
      //x += fTPCCenter;

      fShotPos = TLorentzVector(fCenter.X(), fCenter.Y(), fCenter.Z(), t);


      // Choose angles
      double costh = flat.fire();
      double sinth = std::sqrt(1.0 - cet::square(costh));
      double phi   = 2 * M_PI * flat.fire();

      // Generate momentum 4-vector

      fShotMom = TLorentzVector( p*sinth*cos(phi),
                                 p*sinth*sin(phi),
                                 p*costh,
                                 p                  );

      int trackid = -(mct.NParticles()+1); // set track id to -i as these are all primary particles and have id <= 0
      std::string primary("primary");
      int PDG=0; //optical photons have PDG 0

      simb::MCParticle part(trackid, PDG, primary);
      part.AddTrajectoryPoint(fShotPos, fShotMom);

      if(fFillTree)
        fPhotonsGenerated->Fill();

      mct.Add(part);
    }
    
    return mct;
  } // LightSource::Sample()

} // namespace evgen

DEFINE_ART_MODULE(evgen::LightSource)
