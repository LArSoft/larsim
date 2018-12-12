////////////////////////////////////////////////////////////////////////
/// \file  MUSUN_module.cc
/// \brief Generator for underground muon propagation
///
/// Module designed to propagate muons underground
///  
/// For a description of how to use the module see DUNE DocDB  
/// It is highly recommended that you read it before use.....
///
/// Original MUSUN code was written by Vitaly A. Kudryavtsev (University of Sheffield).
/// Conversion from Fortran to C was done by Kareem Kazkaz (LLNL) with help from David Woodward (Sheffield)
/// Default slant depths and surface profile work done Martin Richardson (Sheffield)
/// Notable contributions to surface map profiling also made by Chao Zhang (USD) and Jeff de Jong (Oxford)
/// Interfaced with LArSoft by Karl Warburton (University of Sheffield) with necessary changes.
/// ------------START OF NECCESSARY CHANGES---------------
/// Changing variables to fcl parameters
/// Restructuring to fit LArSoft design 
/// Co-ordinate transformation... y -> x, x -> z, z -> y ????SHOULD BE x -> -z AS PER VITALY???? 
/// Restore cavern angle convention to Vitaly's definition.
/// ------------END OF NECCESSARY CHANGES-----------------
/// \author  k.warburton@sheffield.ac.uk
///
/// Further Notes - Taken from Kareem Kazkaz; 
///
/// This C++ code is a port of the musun-surf.f and test-musun-surf.f code written
/// by Vitaly Kudryavtsev of the University of Sheffield. It generates muons with
/// energy, position, and direction specific to the Davis cavern at the Sanford
/// Underground Research Facility in Lead, South Dakota.
///
/// This C++ code was ported by Kareem Kazkaz, kareem@llnl.gov, (925) 422-7208
///
/// Here are the notes from Vitaly:
///
///c	The code samples single atmospheric muons at the SURF
///c       underground laboratory (Davis' cavern)
///c       (taking into account the slant depth distribution)
///c	in the energy range E1-E2, zenith angle range theta1-theta2 (0-90 degrees)
///c	and azimuthal angle range phi1-phi2 (0-360 degrees).
///c       At present only the following ranges of parameters are supported:
///c       E1 = 1 GeV, E2 = 10^6 GeV, theta1 = 0, theta2 = 90 deg, phi1 = 0, phi2 = 360 deg.
///c
///c	Program uses muon energy spectra at various depths and zenith
///c	angles obtained with MUSIC code for muon propagation and Gaisser's
///c	formula for muon spectrum at sea level
///c	(T.K.Gaisser, Cosmic Rays and Particle Physics, Cambridge
///c	University Press, 1990) modified for large zenith angles and
///c	prompt muon flux with normalisation and spectral index
///c       that fit LVD data: gamma = 2.77, A = 0.14.
///c       Density of rock is assumed to be 2.70 g/cm^3 but can be changed
///c       during initialisation (previous step, ask the author).
///c
///c       Muon flux through a sphere (Chao's coordinates) = 6.33x10^(-9) cm^(-2) s^(-1) (gamma=2.77) - old
///c       Muon flux through a sphere (Martin's coordinates) = 6.16x10^(-9) cm^(-2) s^(-1) (gamma=2.77) - new
///c       Muon flux through the cuboid (30x22x24 m^3) = 0.0588 s^(-1) (gamma=2.77)
///c
///c	Note: the muon spectrum at sea level does not take into account
///c	the change of the primary spectrum slope at and above the knee
///c	region (3*10^15-10^16 eV).
///c
///c	Program uses the tables of muon energy spectra at various
///c       zenith and azimuthal angles at SURF
///c       calculated with the muon propagation code MUSIC and the
///c       angular distribution of muon intensities at SURF (4850 ft level).
///c
///c       Coordinate system for the muon intensities
///c       is determined by the mountain profile provided
///c       by Chao Zhang (USD, South Dakota): x-axis is pointing to the East.
///c       Muons are sampled on a surface of a rectangular parallelepiped,  
///c       the planes of which are parallel to the walls of the cavern.
///c       The long side of the cavern is pointing at 6.6 deg from the North
///c       to the East (or 90 - 6.6 deg from the East to the North).
///c       Muon coordinates and direction cosines are then given in the
///c       coordinate system related to the cavern with x-axis
///c       pointing along the long side of the cavern at 6.6 deg from the
///c       North to the East.
///c       The angle phi is measured from the positive x-axis towards
///c       positive y-axis (0-360 deg).
///c       Z-axis is pointing upwards. 
///
/// Further notes taken from Vitaly's original code, as written for LBNE 
/// (note the differing definition for theta, which is used here,
/// theta is measured from East to South, not North to East.
///c
///c	The code samples single atmospheric muons at the LBNE
///c       underground site (taking into account the slant depth
///c       distribution calculated by Martin Richardson, Sheffield)
///c	in the energy range E1-E2, zenith angle range theta1-theta2 (0-90 degrees)
///c	and azimuthal angle range phi1-phi2 (0-360 degrees).
///c       At present only the following ranges of parameters are supported:
///c       E1 = 1 GeV, E2 = 10^6 GeV, theta1 = 0, theta2 = 90 deg, phi1 = 0, phi2 = 360 deg.
///c
///c       The LBNE far detector site has coordinates:
///c       Latitude    = 44° 20' 45.21" N
///c       Longitude = 103° 45' 16.13" W
///c       Elevation = 355.8 ft (108.4 m)
///c       from e-mail from Virgil Bocean on 30 May 2013
///c       confirmed by Tracy Lundin on 16 December 2013
///c
///c       For slant depth distribution and muon intensities the azimuthal
///c       angle phi is calculated from East to North.
///c       The long side of the cavern is assumed to be pointing to the
///c       beam (Fermilab) at an angle of 7 deg to the South from the
///c       East.
///c       Muon direction cosines are given in the coordinate system
///c       related to the cavern with positive x-axis pointing to the beam
///c       so direction cosine wrt x-axis is -1 if a muon is coming in the
///c       same direction (wrt to x-axis) as the neutrino beam.
///c
///c
///c	Program uses muon energy spectra at various depths and zenith
///c	angles obtained with MUSIC code for muon propagation and Gaisser's
///c	formula for muon spectrum at sea level
///c	(T.K.Gaisser, Cosmic Rays and Particle Physics, Cambridge
///c	University Press, 1990) modified for large zenith angles and
///c	prompt muon flux with parameters from the best fit to the LVD data. 
///c
///c       Muon flux through a sphere = 6.65x10^(-9) cm^(-2) s^(-1) (gamma=2.77) - old
///c       Muon flux through a sphere = 6.65x10^(-9) cm^(-2) s^(-1) (gamma=2.77) - new
///c       Muon flux through the cuboid (100x40x50 m^3) = 0.3074 s^(-1) (gamma=2.77)
///c
///c	Note: the muon spectrum at sea level does not take into account
///c	the change of the primary spectrum slope at and above the knee
///c	region (3*10^15-10^16 eV).
///c
///c	Program uses the tables of muon energy spectra at various
///c       zenith and azimuthal angles at DUSEL
///c       calculated with the muon propagation code MUSIC and the
///c       angular distribution of muon intensities at 4850 ft level
///c       (location of the LBNE far detector).
///c
///c       Coordinate system is determined by the mountain profile provided
///c       by Jeff de Jong (Cambridge) and Chao Zhang (USD, South Dakota).
///c       This is done to allow the simulation of muons on a surface of
///c       a rectangular parallelepiped, the planes of which are parallel
///c       to the walls of the laboratory halls.
///c       The angle phi in the slant depth distribution file is measured
///c       from the positive x-axis towards positive y-axis.
///c       Z-axis is pointing upwards. 
///c       Note that I use here the azimuthal angle range from 0 to 360 deg.
///c
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_MUSUN
#define EVGEN_MUSUN

// C++ includes.
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <memory>
#include <iterator>
#include <vector>
#include <utility> 
#include <sys/stat.h>
#include <exception>
#include <map>
#include <vector>
#include <algorithm>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"

#include "TTree.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TF2.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TAxis.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// for c2: NTPCs is no longer used
//const int NTPCs = 300;


namespace simb { class MCTruth; }

namespace evgen {

  /// module to produce single or multiple specified particles in the detector
  class MUSUN : public art::EDProducer {

  public:
    explicit MUSUN(fhicl::ParameterSet const& pset);
    virtual ~MUSUN();

    // This is called for each event.
    void produce(art::Event& evt);
    void beginJob();
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);
    void endRun(art::Run& run);

   private:

    void SampleOne(unsigned int   i, simb::MCTruth &mct, CLHEP::HepRandomEngine& engine);
 
    void initialization( double theta1, double theta2, double phi1, double phi2,
			 int figflag, double s_hor, double s_ver1, double s_ver2, double &FI );

    void sampling( double &E, double &theta, double &phi, double &dep, CLHEP::HepRandomEngine& engine );


    static const int kGAUS = 1;    

    int                 fPDG;            ///< PDG code of particles to generate   
    double              fChargeRatio;    ///< Charge ratio of particle / anti-particle  

    std::string         fInputDir;       ///< Input Directory
    std::string         fInputFile1;     ///< Input File 1
    std::string         fInputFile2;     ///< Input File 2
    std::string         fInputFile3;     ///< Input File 3

    double              fCavernAngle;    ///< Angle of the detector from the North to the East.
    double              fRockDensity;    ///< Default rock density is 2.70 g cm-3. If this is
                                         ///< changed then the three input files need to be 
                                         ///< remade. If there is a desire for this contact
                                         ///< Vitaly Kudryavtsev at V.Kudryavtsev@shef.ac.uk

    double              fEmin;           ///< Minimum Kinetic Energy (GeV)
    double              fEmax;           ///< Maximum Kinetic Energy (GeV)

    double              fThetamin;       ///< Minimum theta
    double              fThetamax;       ///< Maximum theta
    double              fPhimin;         ///< Minimum phi
    double              fPhimax;         ///< Maximum phi

    int                 figflag;          ///< If want sampled from sphere or parallelepiped
    double              fXmin;           ///< Minimum X position
    double              fYmin;           ///< Minimum Y position
    double              fZmin;           ///< Minimum Z position
    double              fXmax;           ///< Maximum X position
    double              fYmax;           ///< Maximum Y position
    double              fZmax;           ///< Maximum Z position

    double              fT0;             ///< Central t position (ns) in world coordinates
    double              fSigmaT;         ///< Variation in t position (ns)  
    int                 fTDist;          ///< How to distribute t  (gaus, or uniform)
       
    //Define TFS histograms.....
    /*
    TH1D* hPDGCode;
    TH1D* hPositionX;
    TH1D* hPositionY;
    TH1D* hPositionZ;
    TH1D* hTime;
    TH1D* hMomentumHigh; 
    TH1D* hMomentum; 
    TH1D* hEnergyHigh;
    TH1D* hEnergy;
    TH1D* hDepth;
    TH1D* hDirCosineX; 
    TH1D* hDirCosineY; 
    TH1D* hDirCosineZ;
    TH1D* hTheta;
    TH1D* hPhi;
    */
    int PdgCode;
    double Energy, phi, theta, dep, Time;
    double Momentum, px0, py0, pz0;
    double x0, y0, z0, cx, cy, cz;

    //Define some variables....
    double FI      = 0.;
    double s_hor   = 0.;
    double s_ver1  = 0.;
    double s_ver2  = 0.;

    double spmu[121][62][51];
    double fnmu[32401];
    double depth[360][91];
    double fmu[360][91];
    // for c2: e1 and e2 are unused
    //double e1, e2, the1, the2, ph1, ph2;
    double the1, the2, ph1, ph2;
    double se = 0.;
    double st = 0.;
    double sp = 0.;
    double sd = 0.;

    unsigned int NEvents = 0;

    // TTree
    TTree* fTree;
    /*
    // TTree for CryoPos
    TTree* fCryos;
    int NumTPCs;
    double TPCMinX[NTPCs];
    double TPCMaxX[NTPCs];
    double TPCMinY[NTPCs];
    double TPCMaxY[NTPCs];
    double TPCMinZ[NTPCs];
    double TPCMaxZ[NTPCs];
    double CryoSize[6];
    double DetHall[6];
    */
  };
}

namespace evgen{

  //____________________________________________________________________________
  MUSUN::MUSUN(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
  {
    this->reconfigure(pset);

    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<rndm::NuRandomService>()
      ->createEngine(*this, pset, "Seed");

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();
  }
  //____________________________________________________________________________
  MUSUN::~MUSUN()
  {
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  //  Reconfigure
  ////////////////////////////////////////////////////////////////////////////////
  void MUSUN::reconfigure(fhicl::ParameterSet const& p)
  {
    // do not put seed in reconfigure because we don't want to reset 
    // the seed midstream

    fPDG           = p.get< int                 >("PDG");
    fChargeRatio   = p.get< double              >("ChargeRatio");

    fInputDir      = p.get< std::string         >("InputDir");
    fInputFile1    = p.get< std::string         >("InputFile1");
    fInputFile2    = p.get< std::string         >("InputFile2");
    fInputFile3    = p.get< std::string         >("InputFile3");

    fCavernAngle   = p.get< double              >("CavernAngle");
    fRockDensity   = p.get< double              >("RockDensity");

    fEmin          = p.get< double              >("Emin");
    fEmax          = p.get< double              >("Emax");

    fThetamin      = p.get< double              >("Thetamin");
    fThetamax      = p.get< double              >("Thetamax");

    fPhimin        = p.get< double              >("Phimin");
    fPhimax        = p.get< double              >("Phimax");

    figflag        = p.get<int                  >("igflag");
    fXmin          = p.get<double               >("Xmin");
    fYmin          = p.get<double               >("Ymin");
    fZmin          = p.get<double               >("Zmin");
    fXmax          = p.get<double               >("Xmax");
    fYmax          = p.get<double               >("Ymax");
    fZmax          = p.get<double               >("Zmax");

    fT0            = p.get< double              >("T0");
    fSigmaT        = p.get< double              >("SigmaT");
    fTDist         = p.get< int                 >("TDist");

    return;
  }
  ////////////////////////////////////////////////////////////////////////////////
  //  Begin Job
  ////////////////////////////////////////////////////////////////////////////////
  void MUSUN::beginJob()
  {
    // Make the Histograms....
    art::ServiceHandle<art::TFileService> tfs;
    /*
    hPDGCode      = tfs->make<TH1D>("hPDGCode"     ,"PDG Code"      ,30 , -15           , 15             );
    hPositionX    = tfs->make<TH1D>("hPositionX"   ,"Position (cm)" ,500, ( fXmin - 10 ), ( fXmax + 10 ) );
    hPositionY    = tfs->make<TH1D>("hPositionY"   ,"Position (cm)" ,500, ( fYmin - 10 ), ( fYmax + 10 ) );
    hPositionZ    = tfs->make<TH1D>("hPositionZ"   ,"Position (cm)" ,500, ( fZmin - 10 ), ( fZmax + 10 ) );
    hTime         = tfs->make<TH1D>("hTime"        ,"Time (s)"      ,500, 0    , 1e6       );
    hMomentumHigh = tfs->make<TH1D>("hMomentumHigh","Momentum (GeV)",500, fEmin, fEmax     ); 
    hMomentum     = tfs->make<TH1D>("hMomentum"    ,"Momentum (GeV)",500, fEmin, fEmin+1e3 );
    hEnergyHigh   = tfs->make<TH1D>("hEnergyHigh"  ,"Energy (GeV)"  ,500, fEmin, fEmax     );
    hEnergy       = tfs->make<TH1D>("hEnergy"      ,"Energy (GeV)"  ,500, fEmin, fEmin+1e3 );
    hDepth        = tfs->make<TH1D>("hDepth"       ,"Depth (m)"     ,800, 0    , 14000     );

    hDirCosineX = tfs->make<TH1D>("hDirCosineX","Normalised Direction cosine",500, -1, 1 ); 
    hDirCosineY = tfs->make<TH1D>("hDirCosineY","Normalised Direction cosine",500, -1, 1 ); 
    hDirCosineZ = tfs->make<TH1D>("hDirCosineZ","Normalised Direction cosine",500, -1, 1 );

    hTheta      = tfs->make<TH1D>("hTheta"     ,"Angle (degrees)",500, 0, 90  );
    hPhi        = tfs->make<TH1D>("hPhi"       ,"Angle (degrees)",500, 0, 365 );
    */
    fTree = tfs->make<TTree>("Generator","analysis tree");
    fTree->Branch("particleID",&PdgCode, "particleID/I");
    fTree->Branch("energy" ,&Energy , "energy/D");
    fTree->Branch("time"   ,&Time   , "Time/D"  );
    fTree->Branch("posX"   ,&x0     , "posX/D"  );
    fTree->Branch("posY"   ,&y0     , "posY/D"  );
    fTree->Branch("posZ"   ,&z0     , "posZ/D"  );
    fTree->Branch("cosX"   ,&cx     , "cosX/D"  );
    fTree->Branch("cosY"   ,&cy     , "cosY/D"  );
    fTree->Branch("cosZ"   ,&cz     , "cosZ/D"  );
    fTree->Branch("theta"  ,&theta  , "theta/D" );
    fTree->Branch("phi"    ,&phi    , "phi/D"   );
    fTree->Branch("depth"  ,&dep    , "dep/D"   );
    /*
    fCryos = tfs->make<TTree>("CryoSizes","cryo tree");
    fCryos->Branch("NumTPCs" , &NumTPCs , "NumTPCs/I"       );
    fCryos->Branch("TPCMinX" , &TPCMinX , "TPCMinX[NumTPCs]/D");
    fCryos->Branch("TPCMaxX" , &TPCMaxX , "TPCMaxX[NumTPCs]/D");
    fCryos->Branch("TPCMinY" , &TPCMinY , "TPCMinY[NumTPCs]/D");
    fCryos->Branch("TPCMaxY" , &TPCMaxY , "TPCMaxY[NumTPCs]/D");
    fCryos->Branch("TPCMinZ" , &TPCMinZ , "TPCMinZ[NumTPCs]/D");
    fCryos->Branch("TPCMaxZ" , &TPCMaxZ , "TPCMaxZ[NumTPCs]/D");
    fCryos->Branch("CryoSize", &CryoSize, "CryoSize[6]/D"   );
    fCryos->Branch("DetHall" , &DetHall , "DetHall[6]/D"    );
    */
 }

  ////////////////////////////////////////////////////////////////////////////////
  //  Begin Run
  ////////////////////////////////////////////////////////////////////////////////
  void MUSUN::beginRun(art::Run& run)
  {
    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
    
    // Check fcl parameters were set correctly
    if ( fThetamax > 90.5 ) throw cet::exception("MUSUNGen") << "\nThetamax has to be less than " << M_PI/2 << ", but was entered as " << fThetamax << ", this causes an error so leaving program now...\n\n";  
    if ( fThetamin < 0    ) throw cet::exception("MUSUNGen") << "\nThetamin has to be more than 0, but was entered as " << fThetamin << ", this causes an error so leaving program now...\n\n";      
    if ( fThetamax < fThetamin ) throw cet::exception("MUSUNGen") << "\nMinimum angle is bigger than maximum angle....causes an error so leaving program now....\n\n"; 
    if ( fPhimax > 360.5  ) throw cet::exception("MUSUNGen") << "\nPhimax has to be less than " << 2*M_PI << ", but was entered as " << fPhimax << ", this cause an error so leaving program now...\n\n"; 
    if ( fPhimin < 0      ) throw cet::exception("MUSUNGen") << "\nPhimin has to be more than 0, but was entered as " << fPhimin << ", this causes an error so leaving program now...\n\n"; 
    if ( fPhimax < fPhimin) throw cet::exception("MUSUNGen") << "\nMinimum angle is bigger than maximum angle....causes an error so leaving program now....\n\n";
    if ( fEmax   < fEmin  ) throw cet::exception("MUSUNGen") << "\nMinimum energy is bigger than maximum energy....causes an error so leaving program now....\n\n";
     
    
    run.put(std::move(runcol));

    // area of the horizontal plane of the parallelepiped
    s_hor = (fZmax-fZmin)*(fXmax-fXmin);
    // area of the vertical plane of the parallelepiped, perpendicular to z-axis
    s_ver1 = (fXmax-fXmin)*(fYmax-fYmin);
    // area of the vertical plane of the parallelepiped, perpendicular to x-axis
    s_ver2 = (fZmax-fZmin)*(fYmax-fYmin);
    
    //std::cout << s_hor << " " << s_ver1 << " " << s_ver2 << std::endl;

    initialization(fThetamin,fThetamax,fPhimin,fPhimax,figflag,s_hor,s_ver1,s_ver2,FI );
    
    std::cout << "Material - SURF rock" << std::endl;
    std::cout << "Density = " << fRockDensity << " g/cm^3" << std::endl;
    std::cout << "Parameters for muon spectrum are from LVD best fit" << std::endl;
    std::cout << "Muon energy range = " << fEmin << " - " << fEmax << " GeV" << std::endl;
    std::cout << "Zenith angle range = " << fThetamin << " - " << fThetamax << " degrees" << std::endl;
    std::cout << "Azimuthal angle range = " << fPhimin << " - " << fPhimax << " degrees" << std::endl;
    std::cout << "Global intensity = " << FI << " (cm^2 s)^(-1) or s^(-1) (for muons on the surface)" << std::endl;
    /*
    NumTPCs = geo->NTPC(0);
    std::cout << "There are " << NumTPCs << " in cryostat 0" << std::endl;
    for (unsigned int c=0; c<geo->Ncryostats(); c++) {
      const geo::CryostatGeo& cryostat=geo->Cryostat(c);
      geo->CryostatBoundaries( CryoSize, 0 );
      std::cout << "Cryo bounds " << CryoSize[0] << " "<< CryoSize[1] << " "<< CryoSize[2] << " "<< CryoSize[3] << " "<< CryoSize[4] << " "<< CryoSize[5] << std::endl;
      for (unsigned int t=0; t<cryostat.NTPC(); t++) {
	geo::TPCID id;
	id.Cryostat=c;
	id.TPC=t;
	id.isValid=true;
	const geo::TPCGeo& tpc=cryostat.TPC(id);
	TPCMinX[t] = tpc.MinX();
	TPCMaxX[t] = tpc.MaxX();
	TPCMinY[t] = tpc.MinY();
	TPCMaxY[t] = tpc.MaxY();
	TPCMinZ[t] = tpc.MinZ();
	TPCMaxZ[t] = tpc.MaxZ();
	std::cout << t << "\t" << TPCMinX[t] << " " << TPCMaxX[t] << " " << TPCMinY[t] << " " << TPCMaxY[t] << " " << TPCMinZ[t] << " " << TPCMaxZ[t] << std::endl;
      }
    }
    fCryos -> Fill();
    */
    return;
  }

  ////////////////////////////////////////////////////////////////////////////////
  //  End Run
  ////////////////////////////////////////////////////////////////////////////////
  void MUSUN::endRun(art::Run& run)
  {
    std::cout << "\n\nNumber of muons = " << NEvents << std::endl;
    std::cout << "Mean muon energy = " << se/NEvents << " GeV" << std::endl;
    std::cout << "Mean zenith angle (theta) = " << st/NEvents << " degrees" << std::endl;
    std::cout << "Mean azimuthal angle (phi)= " << sp/NEvents << " degrees" << std::endl;
    std::cout << "Mean slant depth = " << sd/NEvents << " m w.e." << std::endl;
  }
  ////////////////////////////////////////////////////////////////////////////////
  //  Produce
  ////////////////////////////////////////////////////////////////////////////////  
  void MUSUN::produce(art::Event& evt)
  {
    ///unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    
    ++NEvents;

    // get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine(art::ScheduleID::first(),
                                                    moduleDescription().moduleLabel());
    SampleOne(NEvents, truth, engine);

    LOG_DEBUG("MUSUN") << truth;

    truthcol->push_back(truth);

    evt.put(std::move(truthcol));

    return;
  }
  ////////////////////////////////////////////////////////////////////////////////
  //  Sample One
  ////////////////////////////////////////////////////////////////////////////////
  // Draw the type, momentum and position of a single particle from the 
  // FCIHL description
  void MUSUN::SampleOne(unsigned int i, simb::MCTruth &mct, CLHEP::HepRandomEngine& engine)
  {
    CLHEP::RandFlat   flat(engine);
    CLHEP::RandGaussQ gauss(engine);

    Energy = 0;
    theta  = 0;
    phi    = 0; 
    dep    = 0;
    Time   = 0;
    
    sampling( Energy, theta, phi, dep, engine );

    theta = theta* M_PI/180;

    //  changing the angle phi so z-axis is positioned along the long side
    //  of the cavern pointing at 14 deg from the North to the East.
    //  phi += (90. - 14.0);
    //  Want our co-ord rotation going from East to South.
    phi += fCavernAngle;
    if( phi >= 360. )
      phi -= 360.;
    if( phi < 0 )
      phi += 360.;
    phi *= M_PI / 180.;
    
    // set track id to -i as these are all primary particles and have id <= 0
    int trackid = -1*(i+1);
    std::string primary("primary");
    
    // Work out whether particle/antiparticle, and mass...
    double m = 0.0;
    PdgCode = fPDG;
    double ChargeCheck = 1./ ( 1 + fChargeRatio );
    double pdgfire = flat.fire();
    if ( pdgfire < ChargeCheck ) PdgCode=-PdgCode;

    static TDatabasePDG  pdgt;
    TParticlePDG* pdgp = pdgt.GetParticle(PdgCode);
    if (pdgp) m = pdgp->Mass();

    //std::cout << pdgfire << " " << ChargeCheck << " " << PdgCode << " " << m << std::endl;
    
    // Work out T0...
    if(fTDist==kGAUS){
      Time = gauss.fire(fT0, fSigmaT);
    }
    else {
      Time = fT0 + fSigmaT*(2.0*flat.fire()-1.0);
    }
	
    //  The minus sign above is for y-axis pointing up, so the y-momentum
    //  is always pointing down
    cx       = -sin(theta)*sin(phi);
    cy       = -cos(theta);
    cz       = +sin(theta)*cos(phi);
    Momentum = std::sqrt(Energy*Energy-m*m); // Get momentum
    px0      = Momentum * cx;
    py0      = Momentum * cy;
    pz0      = Momentum * cz;
    TLorentzVector pvec(px0, py0, pz0, Energy );   
    
    //  Muon coordinates
    double sh1 = s_hor * cos(theta);
    double sv1 = s_ver1 * sin(theta) * fabs(cos(phi));
    double sv2 = s_ver2 * sin(theta) * fabs(sin(phi));
    double ss = sh1 + sv1 + sv2;
    double xfl1 = flat.fire();
    if( xfl1 <= sh1/ss ) {
      x0 = (fXmax - fXmin)*flat.fire() + fXmin;
      y0 = fYmax;
      z0 = (fZmax - fZmin)*flat.fire() + fZmin;
    } else if( xfl1 <= (sh1+sv1)/ss ) {
      x0 = (fXmax - fXmin)*flat.fire() + fXmin;
      y0 = (fYmax - fYmin)*flat.fire() + fYmin;
      if( cz >= 0 ) z0 = fZmax;
      else z0 = fZmin;
    } else {
      if( cx >= 0 ) x0 = fXmin;
      else x0 = fXmax;
      y0 = (fYmax - fYmin)*flat.fire() + fYmin;
      z0 = (fZmax - fZmin)*flat.fire() + fZmin;
    }
    // Make Lorentz vector for x and t....
    TLorentzVector pos(x0, y0, z0, Time); 
   
    //  Parameters written to the file muons_surf_v2_test*.dat
    //      nmu - muon sequential number
    //      id_part - muon charge (10 - positive, 11 - negative )
    //      Energy  - total muon energy in GeV assuming ultrarelativistic muons
    //      x0, y0, z0 - muon coordinates on the surface of parallelepiped
    //          specified above; x-axis and y-axis are pointing in the
    //          directions such that the angle phi (from the slant depth
    //          distribution files) is measured from x to y. z-axis is
    //          pointing upwards.
    //      cx, cy, cz - direction cosines.
  
    simb::MCParticle part(trackid, PdgCode, primary);
    part.AddTrajectoryPoint(pos, pvec);

    mct.Add(part);

    theta = theta * 180/M_PI;
    phi   = phi   * 180/M_PI;

    // Sum energies, angles, depth for average outputs. 
    se += Energy;
    st += theta;
    sp += phi;
    sd += dep;
    
    // Fill Histograms.....
    /*
    hPDGCode      ->Fill (PdgCode);
    hPositionX    ->Fill (x0);
    hPositionY    ->Fill (y0);
    hPositionZ    ->Fill (z0);
    hTime         ->Fill (Time);
    hMomentumHigh ->Fill (Momentum); 
    hMomentum     ->Fill (Momentum); 
    hEnergyHigh   ->Fill (Energy);
    hEnergy       ->Fill (Energy);
    hDepth        ->Fill (dep);
    hDirCosineX   ->Fill (cx); 
    hDirCosineY   ->Fill (cy); 
    hDirCosineZ   ->Fill (cz);
    hTheta        ->Fill (theta);
    hPhi          ->Fill (phi);
    */
    /*
    // Write event by event outsputs.....
    std::cout << "Theta: " << theta << " Phi: " << phi << " Energy: " << Energy << " Momentum: " << Momentum << std::endl; 
    std::cout << "x: " <<  pos.X() << " y: " << pos.Y() << " z: " << pos.Z() << " time: " << pos.T() << std::endl;   
    std::cout << "Px: " <<  pvec.Px() << " Py: " << pvec.Py() << " Pz: " << pvec.Pz() << std::endl;
    std::cout << "Normalised..." << cx << " " << cy << " " << cz << std::endl;
    */
    fTree->Fill();

  }

  ////////////////////////////////////////////////////////////////////////////////
  //  initialization
  ////////////////////////////////////////////////////////////////////////////////
  void MUSUN::initialization( double theta1, double theta2, double phi1, double phi2,
		       int figflag, double s_hor, double s_ver1, double s_ver2, double &FI )
  {
    //
    //  Read in the data files
    //
    int lineNumber = 0, index = 0;
    char inputLine[10000];
    std::string fROOTfile; 

    for (int a=0;a<121;++a) for (int b=0;b<62;++b) for (int c=0;c<50;++c) spmu[a][b][c]=0;
    for (int a=0;a<23401;++a)                      fnmu[a]    = 0;
    for (int a=0;a<360;++a) for (int b=0;b<91;++b) depth[a][b]= 0;
    for (int a=0;a<360;++a) for (int b=0;b<91;++b) fmu[a][b]  = 0;

    std::ostringstream File1LocStream;
    File1LocStream  << fInputDir << fInputFile1;
    std::string File1Loc = File1LocStream. str();
    cet::search_path sp1("FW_SEARCH_PATH");
    if( sp1.find_file(fInputFile1, fROOTfile) ) File1Loc = fROOTfile;
    std::ifstream file1( File1Loc.c_str(), std::ios::in );
    if (!file1.good() ) throw cet::exception("MUSUNGen") << "\nFile1 " << fInputFile1 << " not found in FW_SEARCH_PATH or at " <<  fInputDir <<"\n\n";
    
    while( file1.good() ) {
      //std::cout << "Looking at file 1...." << std::endl;
      file1.getline( inputLine, 9999 );
      char *token;
      token = strtok( inputLine, " " );
      while( token != NULL ) {
	//std::cout << "While loop file 1..." << std::endl;
	fmu[index][lineNumber] = atof( token );
	token = strtok( NULL, " " );
	index++;
	if( index == 360 ) {
	  //std::cout << "If statement file 1..." << std::endl;
	  index = 0;
	  lineNumber++;
	}
      }
    }
    file1.close();
    
    std::ostringstream File2LocStream;
    File2LocStream  << fInputDir << fInputFile2;
    std::string File2Loc = File2LocStream. str();
    cet::search_path sp2("FW_SEARCH_PATH");
    if( sp2.find_file(fInputFile2, fROOTfile) ) File2Loc = fROOTfile;
    std::ifstream file2( File2Loc.c_str(), std::ios::binary|std::ios::in );
    if (!file2.good() ) throw cet::exception("MUSUNGen") << "\nFile2 " << fInputFile2 << " not found in FW_SEARCH_PATH or at " <<  fInputDir <<"\n\n";
  
    int i1 = 0, i2 = 0, i3 = 0;
    float readVal;
    while( file2.good() ) {
      //std::cout << "Looking at file 2...." << std::endl;
      file2.read((char *)(&readVal), sizeof(float));
      spmu[i1][i2][i3] = readVal;
      i1++;
      if( i1 == 121 ) {
	//std::cout << "First if statement file 2..." << std::endl;
	i2++;
	i1 = 0;
      }
      if( i2 == 62 ) {
	//std::cout << "Second if statement file 2..." << std::endl;
	i3++;
	i2 = 0;
      }
    }
    file2.close();
    for( int i=0; i<120; i++ )
      for( int j=0; j<62; j++ )
	for( int k=0; k<51; k++ )
	  spmu[i][j][k] = spmu[i+1][j][k];
    spmu[1][1][0] = 0.000853544;
    //std::cout << "Set spmu to some value..." << std::endl;
    
    std::ostringstream File3LocStream;
    File3LocStream  << fInputDir << fInputFile3;
    std::string File3Loc = File3LocStream. str();
    cet::search_path sp3("FW_SEARCH_PATH");
    if( sp3.find_file(fInputFile3, fROOTfile) ) File3Loc = fROOTfile;
    std::ifstream file3( File3Loc.c_str(), std::ios::in );
    if (!file3.good() ) throw cet::exception("MUSUNGen") << "\nFile3 " << fInputFile3 << " not found in FW_SEARCH_PATH or at " <<  fInputDir <<"\n\n";

    lineNumber = index = 0;
    while( file3.good() ) {
      //std::cout << "Looking at file 3...." << std::endl;
      file3.getline( inputLine, 9999 );
      char *token;
      token = strtok( inputLine, " " );
      while( token != NULL ) {
	//std::cout << "While loop file 3..." << std::endl;
	depth[index][lineNumber] = atof( token );
	token = strtok( NULL, " " );
	index++;
	if( index == 360 ) {
	  //std::cout << "If statement file 3..." << std::endl;
	  index = 0;
	  lineNumber++;
	}
      }
    }
    file3.close();
    
    //
    //  Set up variables
    //
    
    the1 = theta1;
    the2 = theta2;
    // for c2: c1 and c2 are unused
    //double c1 = cos(M_PI/180.*theta1);
    //double c2 = cos(M_PI/180.*theta2);
    ph1 = M_PI/180.*phi1;
    ph2 = M_PI/180.*phi2;
    // for c2: dph is unused
    //double dph = ph2-ph1;
    
    int ipc = 1;
    double theta = theta1;
    double dc = 1.;
    double sc = 0.;
    int iteration = 0;
    while( theta < theta2-dc/2. ) {
        theta += dc/2.;
        double theta0 = M_PI/180. * theta;
        double cc = cos(theta0);
        double ash = s_hor * cc;
        double asv01 = s_ver1 * sqrt(1. - cc*cc);
        double asv02 = s_ver2 * sqrt(1. - cc*cc);
        int ic1 = (theta + 0.999);
        int ic2 = ic1 + 1;
        if( ic2 > 91 ) ic2 = 91;
        if( ic1 < 1 ) ic1 = 1;
        double phi = phi1;
        double dp = 1.;
        
        while( phi < phi2-dp/2. ) {
            phi += dp/2.;
            //  the long side of the cavern is pointing at 14 deg to the north:
            //  double phi0 = M_PI / 180. * (phi + 90. - 14);
	    
	    //  Want our co-ord system going from East to South.
	    double phi0 = M_PI / 180. * (phi + fCavernAngle);
            
	    double asv1 = asv01 * fabs(cos(phi0));
            double asv2 = asv02 * fabs(sin(phi0));
            double asv0 = ash + asv1 + asv2;
            double fl = 1.;
            if( figflag == 1 )
                fl = asv0;
            int ip1 = (phi + 0.999);
            int ip2 = ip1 + 1;
            if( ip2 > 360 ) ip2 = 1;
            if( ip1 < 1 ) ip1 = 360;
            double sp1 = 0.;
            
            for( int ii=0; ii<4; ii++ ) {
                int iic = ii/2;
                int iip = ii%2;
		if(ip1==360 && (ii==1 || ii==3) ) iip = -359;
		if( fmu[ip1+iip-1][ic1+iic-1] < 0 ) {
		  if ( pow(10.,fmu[ip1+iip-1][ic1+iic-1]) / 4 > 1e-6 ) {
		    std::cout << "Looking at fmu [ " << ip1 << " + " << iip << " - 1 (" << ip1+iip-1 << ") ] [ " << ic1 << " + " << iic << " - 1 ("<< ic1+iic-1 << ") ] ."
			      << "\nChanging sp1 from " << sp1 << " to " <<  sp1 + pow(10.,fmu[ip1+iip-1][ic1+iic-1]) / 4 << "..........." << sp1 << " + 10 ^ (" << fmu[ip1+iip-1][ic1+iic-1] << ") / 4 "
			      << std::endl;
		  }
		  sp1 = sp1 + pow(10.,fmu[ip1+iip-1][ic1+iic-1]) / 4;
		}
            }
	    /*
	      std::cout << iteration<< " time of new sc value! Theta " << theta << ", phi " << phi + dp / 2. << ", sc = " << sc + sp1 * fl * dp * M_PI / 180. * sin(theta0) * dc * M_PI / 180. << " = " 
	      << sc << " + " << sp1 << " * " << fl << " * " << dp << " * " << M_PI/180 << " * sin(" << theta0 << ") * " << dc << " * " << M_PI/180 << ".....sin(theta)=" << sin(theta) << "\n" 
	      << std::endl; */
	    sc = sc + sp1 * fl * dp * M_PI / 180. * sin(theta0) * dc * M_PI / 180.;
	    ++iteration;
	    ipc = ipc + 1;
            fnmu[ipc-1] = sc;
            phi = phi + dp / 2.;
        }
        
        theta = theta + dc / 2.;
    }
    //std::cout << *FI << " = " << sc << std::endl;
    FI = sc;
    for( int ipc1 = 0; ipc1 < ipc; ipc1++ )
      fnmu[ipc1] = fnmu[ipc1] / fnmu[ipc-1];
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////
  //  sampling
  ////////////////////////////////////////////////////////////////////////////////
  void MUSUN::sampling( double &E, double &theta, double &phi, double &dep, CLHEP::HepRandomEngine& engine )
  {
    CLHEP::RandFlat   flat(engine);
    CLHEP::RandGaussQ gauss(engine);

    #if 0 // this code is disabled for good
    double xfl = flat.fire();
    int loIndex = 0, hiIndex = 32400;
    int i = (loIndex+hiIndex)/2;
    bool foundIndex = false;
    if( xfl < fnmu[loIndex] ) {
      i = loIndex;
      foundIndex = true;
    } else if ( xfl > fnmu[hiIndex] ) {
      i = hiIndex;
      foundIndex = true;
    } else if ( xfl > fnmu[i-1] && xfl <= fnmu[i] )
      foundIndex = true;
    while( !foundIndex ) {
      if( xfl < fnmu[i] )
	hiIndex = i;
      else
	loIndex = i;
      i = (loIndex + hiIndex)/2;
      
      if( xfl > fnmu[i-1] && xfl <= fnmu[i] )
	foundIndex = true;
    }
    #else
    double xfl = flat.fire();
    int i = 0;
    while ( xfl > fnmu[i] ) ++i;
    #endif
    int ic = (i-2)/360;
    int ip = i-2-ic*360;
        
    xfl = flat.fire();
    theta = the1 + ((double)ic+xfl);
    xfl = flat.fire();
    phi = ph1 + ((double)ip+xfl);
    if ( phi > 360 ) phi = phi -360;
    dep = depth[ip][ic] * fRockDensity;
    
    int ic1 = cos(M_PI/180.*theta) * 50.;
    if( ic1 < 0 )
      ic1 = 0;
    if( ic1 > 50 )
      ic1 = 50;
    int ip1 = dep / 200. - 16;
    if( ip1 < 0 )
      ip1 = 0;
    if( ip1 > 61 )
      ip1 = 61;

    xfl = flat.fire();
    #if 0
    loIndex = 0, hiIndex = 120;
    int j = (loIndex+hiIndex)/2;
    foundIndex = false;
    if( xfl < spmu[loIndex][ip1][ic1] ) {
      j = loIndex;
      foundIndex = true;
    } else if ( xfl > spmu[hiIndex][ip1][ic1] ) {
      j = hiIndex;
      foundIndex = true;
    } else if ( xfl > spmu[j-1][ip1][ic1] && xfl <= spmu[j][ip1][ic1] )
      foundIndex = true;
    while( !foundIndex ) {
      if( xfl < spmu[j][ip1][ic1] )
	hiIndex = j;
      else
	loIndex = j;
      j = (loIndex + hiIndex)/2;
      
      if( xfl > spmu[j-1][ip1][ic1] && xfl <= spmu[j][ip1][ic1] )
	foundIndex = true;
    }
    #else
    int j = 0;
    while ( xfl > spmu[j][ip1][ic1] ) ++j;
    #endif

    double En1 = 0.05 * (j-1);
    double En2 = 0.05 * (j);
    E = pow(10.,En1 + (En2 - En1)*flat.fire());
    
    return;
  }
  
}//end namespace evgen

namespace evgen{
  
  DEFINE_ART_MODULE(MUSUN)
  
}//end namespace evgen

#endif
////////////////////////////////////////////////////////////////////////
