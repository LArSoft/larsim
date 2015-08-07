////////////////////////////////////////////////////////////////////////
/// \file  SingleGen_plugin.cc
/// \brief Generator for cosmic-rays
///
/// Module designed to produce muons for a MC event using a Gaissers
/// parametisation. 
/// For a description of how to use the module see DUNE DocDB 10741 
/// It is highly reccommended that you read it before use.....
///
/// \version $Id: GaisserParam.cxx,v 1.4 2015/04/20 09:54:01 brebel Exp $
/// \author  k.warburton@sheffield.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_GAISSERPARAM
#define EVGEN_GAISSERPARAM

// C++ includes.
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
#include "cetlib/exception.h"

// art extensions
#include "artextensions/SeedService/SeedService.hh"

// nutools includes
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "EventGeneratorBase/evgenbase.h"

// lar includes
#include "Geometry/Geometry.h"
#include "SummaryData/RunData.h"

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


namespace simb { class MCTruth; }

namespace evgen {

  /// module to produce single or multiple specified particles in the detector
  class GaisserParam : public art::EDProducer {

  public:
    explicit GaisserParam(fhicl::ParameterSet const& pset);
    virtual ~GaisserParam();

    // This is called for each event.
    void produce(art::Event& evt);
    void beginJob();
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);

    // Defining public maps.......
    typedef std::map<double, TH1*> dhist_Map;
    typedef std::map<double, TH1*>::iterator dhist_Map_it;
    TFile* m_File;
    dhist_Map* m_PDFmap;
    TH1* m_thetaHist;  

  private:

    void SampleOne(unsigned int   i, simb::MCTruth &mct);   
    void Sample(simb::MCTruth &mct);        

    std::pair<double,double> GetThetaAndEnergy(double rand1, double rand2);
    void MakePDF();
    void ResetMap();
    double GaisserMuonFlux_Integrand(Double_t *x, Double_t *par);    
    double GaisserFlux(double e, double theta);
    std::vector<double> GetBinning(const TAxis* axis, bool finalEdge=true);


    static const int kGAUS = 1;    

    int                 fMode;           ///< Particle Selection Mode 
                                         ///< 0--generate a list of all particles, 
                                         ///< 1--generate a single particle selected randomly from the list
    bool                fPadOutVectors;  ///< Select to pad out configuration vectors if they are not of 
					 ///< of the same length as PDG  
                                         ///< false: don't pad out - all values need to specified
                                         ///< true: pad out - default values assumed and printed out
    std::vector<int>    fPDG;            ///< PDG code of particles to generate   
    int                 fCharge;         ///< Charge  
    std::string         fInputDir;       ///< Input Directory

    double              fEmin;           ///< Minimum Kinetic Energy (GeV)
    double              fEmax;           ///< Maximum Kinetic Energy (GeV)
    double              fEmid;           ///< Energy to go from low to high (GeV)
    int                 fEBinsLow;       ///< Number of low energy Bins
    int                 fEBinsHigh;      ///< Number of high energy Bins

    double              fThetamin;       ///< Minimum theta
    double              fThetamax;       ///< Maximum theta
    int                 fThetaBins;      ///< Number of theta Bins

    double              fXHalfRange;           ///< Max X position
    double              fYInput;           ///< Max Y position
    double              fZHalfRange;           ///< Max Z position
 
    double              fT0;             ///< Central t position (ns) in world coordinates
    double              fSigmaT;         ///< Variation in t position (ns)  
    int                 fTDist;          ///< How to distribute t  (gaus, or uniform)

    bool                fSetParam;       ///< Which version of Gaissers Param
    bool                fSetRead;        ///< Whether to Read
    bool                fSetWrite;       ///< Whether to Write
    bool                fSetReWrite;     ///< Whether to ReWrite pdfs
    double              fEpsilon;        ///< Minimum integration sum....
       
    //Define TFS histograms.....
    TH1D* fPositionX;
    TH1D* fPositionY;
    TH1D* fPositionZ;
    TH1D* fTime;
    TH1D* fMomentumHigh; 
    TH1D* fMomentum; 
    TH1D* fEnergy;
    TH1D* fDirCosineX; 
    TH1D* fDirCosineY; 
    TH1D* fDirCosineZ;
    TH1D* fTheta;
    TH1D* fPhi;

    //Define some variables....
    double fCryoBoundaries[6];
    double xNeg    = 0;
    double xPos    = 0;
    double zNeg    = 0;
    double zPos    = 0;
    double fCenterX= 0;
    double fCenterZ= 0;
  };
}

namespace evgen{

  //____________________________________________________________________________
  GaisserParam::GaisserParam(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    // create a default random engine; obtain the random seed from SeedService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<artext::SeedService>()
      ->createEngine(*this, pset, "Seed");

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();
  }

  //____________________________________________________________________________
  GaisserParam::~GaisserParam()
  {
  }

  //____________________________________________________________________________
  void GaisserParam::reconfigure(fhicl::ParameterSet const& p)
  {
    // do not put seed in reconfigure because we don't want to reset 
    // the seed midstream

    fPadOutVectors = p.get< bool                >("PadOutVectors");
    fMode          = p.get< int                 >("ParticleSelectionMode");
    fPDG           = p.get< std::vector<int>    >("PDG");

    fCharge        = p.get< int >("Charge");
    fInputDir      = p.get< std::string         >("InputDir");

    fEmin          = p.get< double              >("Emin");
    fEmax          = p.get< double              >("Emax");
    fEmid          = p.get< double              >("Emid");
    fEBinsLow      = p.get< int                 >("EBinsLow");
    fEBinsHigh     = p.get< int                 >("EBinsHigh");

    fThetamin      = p.get< double              >("Thetamin");
    fThetamax      = p.get< double              >("Thetamax");
    fThetaBins     = p.get< int                 >("ThetaBins");
   
    fXHalfRange    = p.get< double              >("XHalfRange");
    fYInput        = p.get< double              >("YInput");
    fZHalfRange    = p.get< double              >("ZHalfRange");

    fT0            = p.get< double              >("T0");
    fSigmaT        = p.get< double              >("SigmaT");
    fTDist         = p.get< int                 >("TDist");

    fSetParam      = p.get< bool                >("SetParam");
    fSetRead       = p.get< bool                >("SetRead");
    fSetWrite      = p.get< bool                >("SetWrite");
    fSetReWrite    = p.get< bool                >("SetReWrite");
    fEpsilon       = p.get< double              >("Epsilon");
    
    return;
  }

  //____________________________________________________________________________
  void GaisserParam::beginJob()
  {
    //Work out center of cryostat(s)
    art::ServiceHandle<geo::Geometry> geom;
    for (unsigned int i=0; i < geom->Ncryostats() ; i++ ) {
      geom->CryostatBoundaries(fCryoBoundaries, i);
      if ( xNeg > fCryoBoundaries[0] ) xNeg = fCryoBoundaries[0];
      if ( xPos < fCryoBoundaries[1] ) xPos = fCryoBoundaries[1];
      if ( zNeg > fCryoBoundaries[4] ) zNeg = fCryoBoundaries[4];
      if ( zPos < fCryoBoundaries[5] ) zPos = fCryoBoundaries[5];
      }
    fCenterX = xNeg + (xPos-xNeg)/2;
    fCenterZ = zNeg + (zPos-zNeg)/2;
    
    // Make the Histograms....
    art::ServiceHandle<art::TFileService> tfs;
    fPositionX    = tfs->make<TH1D>("fPositionX"   ,"Position (cm)" ,500,fCenterX-(fXHalfRange+10) ,fCenterX+(fXHalfRange+10));
    fPositionY    = tfs->make<TH1D>("fPositionY"   ,"Position (cm)" ,500,-(fYInput+10),(fYInput+10));
    fPositionZ    = tfs->make<TH1D>("fPositionZ"   ,"Position (cm)" ,500,fCenterZ-(fZHalfRange+10) ,fCenterZ+(fZHalfRange+10));
    fTime         = tfs->make<TH1D>("fTime"        ,"Time (s)"      ,500,0,1e6);
    fMomentumHigh = tfs->make<TH1D>("fMomentumHigh","Momentum (GeV)",500,0,fEmax); 
    fMomentum     = tfs->make<TH1D>("fMomentum"    ,"Momentum (GeV)",500,0,100);
    fEnergy       = tfs->make<TH1D>("fEnergy"      ,"Energy (GeV)"  ,500,0,fEmax);

    fDirCosineX = tfs->make<TH1D>("fDirCosineX","Normalised Direction cosine",500,-1,1); 
    fDirCosineY = tfs->make<TH1D>("fDirCosineY","Normalised Direction cosine",500,-1,1); 
    fDirCosineZ = tfs->make<TH1D>("fDirCosineZ","Normalised Direction cosine",500,-1,1);

    fTheta      = tfs->make<TH1D>("fTheta"     ,"Angle (radians)",500,-365,365);
    fPhi        = tfs->make<TH1D>("fPhi"       ,"Angle (radians)",500,-365,365);

  }

  //____________________________________________________________________________
  void GaisserParam::beginRun(art::Run& run)
  {
    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));

    // Check fcl parameters were set correctly
    if ( fThetamax > M_PI/2 + 0.01 ) {
      std::cout << "\n\nThetamax has to be less than " << M_PI/2 << ", but was entered as " << fThetamax << ", this cause an error so leaving program now...\n\n" << std::endl; 
      return;
    }
    if ( fThetamin < 0 ) {
      std::cout << "\n\nThetamin has to be more than 0, but was entered as " << fThetamin << ", this cause an error so leaving program now...\n\n" << std::endl; 
      return;
    }
    if ( fThetamax < fThetamin ) {
      std::cout << "\n\nMinimum angle is bigger than maximum angle....causes an error so leaving program now....\n\n" << std::endl;
      return;
    }
    run.put(std::move(runcol));
    MakePDF ();
    
    return;
  }
  
  //____________________________________________________________________________
  void GaisserParam::produce(art::Event& evt)
  {
    ///unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    
    Sample(truth);

    LOG_DEBUG("GaisserParam") << truth;

    truthcol->push_back(truth);

    evt.put(std::move(truthcol));

    return;
  }

  //____________________________________________________________________________
  // Draw the type, momentum and position of a single particle from the 
  // FCIHL description
  void GaisserParam::SampleOne(unsigned int i, simb::MCTruth &mct){

    // get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat   flat(engine);
    CLHEP::RandGaussQ gauss(engine);

    double Time          = 0.0; 
    double Momentum      = 0.0;
    double KineticEnergy = 0.0;
    double Gamma         = 0.0;
    double Energy        = 0.0;
    double DirCosineX    = 0.0;
    double DirCosineY    = 0.0;
    double DirCosineZ    = 0.0;
    double Theta         = 0.0;
    double Phi           = 0.0;
    double pnorm  = 0.0;
    
    TVector3 x;
    TVector3 pmom;

    // set track id to -i as these are all primary particles and have id <= 0
    int trackid = -1*(i+1);
    std::string primary("primary");

    // Work out whether particle/antiparticle, and mass...
    double m =0.0;
    fPDG[i] = 13;
    if (fCharge == 0 ) {
      if(1.0-2.0*flat.fire() > 0) fPDG[i]=-fPDG[i];
    } 
    else if ( fCharge == 1 ) fPDG[i] = 13;
    else if ( fCharge == 2 ) fPDG[i] = -13;
    static TDatabasePDG  pdgt;
    TParticlePDG* pdgp = pdgt.GetParticle(fPDG[i]);
    if (pdgp) m = pdgp->Mass();
    
    // Work out T0...
    if(fTDist==kGAUS){
      Time = gauss.fire(fT0, fSigmaT);
    }
    else {
      Time = fT0 + fSigmaT*(2.0*flat.fire()-1.0);
    }
    
    // Work out Positioning.....
    x[0] = (1.0 - 2.0*flat.fire())*fXHalfRange + fCenterX;
    x[1] = fYInput;
    x[2] = (1.0 - 2.0*flat.fire())*fZHalfRange + fCenterZ;

    // Make Lorentz vector for x and t....
    TLorentzVector pos(x[0], x[1], x[2], Time); 
    
    // Access the pdf map which has been loaded.....
    if(m_PDFmap) {

      //---- get the muon theta and energy from histograms using 2 random numbers      
      std::pair<double,double> theta_energy; //---- muon theta and energy            
      theta_energy = GetThetaAndEnergy(flat.fire(),flat.fire());
      
      //---- Set theta, phi
      Theta = theta_energy.first;  // Angle returned by GetThetaAndEnergy between 0 and pi/2 
      Phi   = M_PI*( 1.0-2.0*flat.fire() ); // Randomly generated angle between -pi and pi

      //---- Set KE, E, p
      KineticEnergy = theta_energy.second; // Energy returned by GetThetaAndEnergy
      Gamma         = 1 + (KineticEnergy/m);
      Energy        = Gamma * m;
      Momentum      = std::sqrt(Energy*Energy-m*m); // Get momentum

      pmom[0] = Momentum*std::sin(Theta)*std::cos(Phi);
      pmom[1] = -Momentum*std::cos(Theta);
      pmom[2] = Momentum*std::sin(Theta)*std::sin(Phi);
      
      pnorm   = std::sqrt( pmom[0]*pmom[0] + pmom[1]*pmom[1] + pmom[2]*pmom[2] );
      DirCosineX = pmom[0] / pnorm;
      DirCosineY = pmom[1] / pnorm;
      DirCosineZ = pmom[2] / pnorm;
    } 
    else {
      std::cout << "MuFlux map hasn't been initialised, aborting...." << std::endl;
      return;
    }    

    TLorentzVector pvec(pmom[0], pmom[1], pmom[2], Energy );

    simb::MCParticle part(trackid, fPDG[i], primary);
    part.AddTrajectoryPoint(pos, pvec);

    fPositionX    ->Fill (x[0]);
    fPositionY    ->Fill (x[1]);
    fPositionZ    ->Fill (x[2]);
    fTime         ->Fill (Time);
    fMomentumHigh ->Fill (Momentum); 
    fMomentum     ->Fill (Momentum); 
    fEnergy       ->Fill (Energy);
    fDirCosineX   ->Fill (DirCosineX); 
    fDirCosineY   ->Fill (DirCosineY); 
    fDirCosineZ   ->Fill (DirCosineZ);
    fTheta        ->Fill (Theta*180/M_PI);
    fPhi          ->Fill (Phi  *180/M_PI);
  
    std::cout << "Theta: " << Theta << " Phi: " << Phi << " KineticEnergy: " << KineticEnergy << " Energy: " << Energy << " Momentum: " << Momentum << std::endl; 
    std::cout << "x: " <<  pos.X() << " y: " << pos.Y() << " z: " << pos.Z() << " time: " << pos.T() << std::endl;   
    std::cout << "Px: " <<  pvec.Px() << " Py: " << pvec.Py() << " Pz: " << pvec.Pz() << std::endl;
    std::cout << "Normalised..." << DirCosineX << " " << DirCosineY << " " << DirCosineZ << std::endl;
            
    mct.Add(part);
  }

  //____________________________________________________________________________
  void GaisserParam::Sample(simb::MCTruth &mct) 
  {
    switch (fMode) {
    case 0: // List generation mode: every event will have one of each
	    // particle species in the fPDG array
      for (unsigned int i=0; i<fPDG.size(); ++i) {
	SampleOne(i,mct);
      }//end loop over particles
      break;
    case 1: // Random selection mode: every event will exactly one particle
            // selected randomly from the fPDG array
      {
	art::ServiceHandle<art::RandomNumberGenerator> rng;
	CLHEP::HepRandomEngine &engine = rng->getEngine();
	CLHEP::RandFlat flat(engine);

	unsigned int i=flat.fireInt(fPDG.size());
	SampleOne(i,mct);
      }
      break;
    default:
      mf::LogWarning("UnrecognizeOption") << "GaisserParam does not recognize ParticleSelectionMode "
					  << fMode;
      break;
    } // switch on fMode

    return;
  }

  //__________________________________
  std::pair<double,double> GaisserParam::GetThetaAndEnergy(double rand1, double rand2)
  {
    if(rand1 < 0 || rand1 > 1) std::cerr << "GetThetaAndEnergy:\tInvalid random number " << rand1 << std::endl;
    if(rand2 < 0 || rand2 > 1) std::cerr << "GetThetaAndEnergy:\tInvalid random number " << rand2 << std::endl;
    
    int thetaBin = 0;
    m_thetaHist->GetBinWithContent(double(rand1),thetaBin,1,m_thetaHist->GetNbinsX(),1.0);  
    if(m_thetaHist->GetBinContent(thetaBin) < rand1) thetaBin += 1;  
    double drand1 = (m_thetaHist->GetBinContent(thetaBin) - rand1)/(m_thetaHist->GetBinContent(thetaBin) - m_thetaHist->GetBinContent(thetaBin-1));  
    double thetaLow = m_thetaHist->GetXaxis()->GetBinLowEdge(thetaBin);
    double thetaUp  = m_thetaHist->GetXaxis()->GetBinUpEdge(thetaBin);  
    double theta = drand1*(thetaUp-thetaLow) + thetaLow;

    int energyBin = 0;
    TH1* energyHist = 0;
    bool notFound = true;
    for(dhist_Map_it mapit=m_PDFmap->begin(); mapit!=m_PDFmap->end(); mapit++){
      if( fabs(mapit->first+thetaLow)<0.000001 ) {
	energyHist = mapit->second;
	notFound = false;
	break;
      }
    }
    if(notFound) std::cout << "GetThetaAndEnergy: ERROR:\tInvalid theta!" << std::endl;
    //  MSG("string = " << To_TString(thetaLow) << ", double = " << thetaLow );

    energyHist->GetBinWithContent(double(rand2),energyBin,1,energyHist->GetNbinsX(),1.0);  
    if(energyHist->GetBinContent(energyBin) < rand2) energyBin += 1;  
    double drand2 = (energyHist->GetBinContent(energyBin) - rand2)/(energyHist->GetBinContent(energyBin) - energyHist->GetBinContent(energyBin-1));  
    double energyLow = energyHist->GetXaxis()->GetBinLowEdge(energyBin);
    double energyUp  = energyHist->GetXaxis()->GetBinUpEdge(energyBin);  
    double energy = drand2*(energyUp-energyLow) + energyLow;
    //  MSG("MuFlux::GetThetaEnergy()\te = " << energy*1000. );

    return std::make_pair(theta,energy);
  }

  //____________________________________________________________________________
  void GaisserParam::MakePDF() 
  {
    std::cout << "In my function MakePDF" << std::endl;
    m_File=0;
    m_PDFmap=0;
    m_thetaHist=0;  
    double TotalMuonFlux=0;

    if(m_PDFmap){
      std::cout << "PDFMAP" << std::endl;
      if(m_PDFmap->size()>0 && !fSetReWrite){
	std::cout << "MakePDF: Map has already been initialised. " << std::endl;
	std::cout << "Do fSetReWrite - true if you really want to override this map." << std::endl;
	return; 
	}
      std::cout << fSetReWrite << std::endl;
      if(fSetReWrite) ResetMap();
    }
    else{
      m_PDFmap = new dhist_Map;
      std::cout << "Making new dhist_Map called m_PDFmap....." << std::endl;
    }
    
    TF2* muonSpec = new TF2("muonSpec", this, 
			    &evgen::GaisserParam::GaisserMuonFlux_Integrand, 
			    fEmin, fEmax, fThetamin, fThetamax, 0, 
			    "GaisserParam", "GaisserMuonFlux_Integrand"
			    );
    //--------------------------------------------
    //------------ Compute the pdfs
    
    //---- compute pdf for the theta
    TotalMuonFlux = muonSpec->Integral(fEmin, fEmax, fThetamin, fThetamax, fEpsilon ); // Work out the muon flux at the surface
    std::cout << "Surface flux of muons = " << TotalMuonFlux << " cm-2 s-1" << std::endl;

    //---- work out if we're reading a file, writing to file, or neither
    std::ostringstream pdfFile;
    pdfFile << "GaisserPDF_"<<fEmin<<"-"<<fEmid<<"-"<<fEmax<<"-"<< fEBinsLow<<"-"<<fEBinsHigh<<"-"<<fThetamin<<"-"<<fThetamax<<"-"<<fThetaBins<<".root"; 
    std::string tmpfileName = pdfFile.str();
    std::replace(tmpfileName.begin(),tmpfileName.end(),'+','0');
    if      (tmpfileName == "GaisserPDF_0-100-100100-1000-10000-0-1.5708-100.root")          tmpfileName = "GaisserPDF_DefaultBins.root";
    else if (tmpfileName == "GaisserPDF_0-100-4000-1000-1000-0-1.5708-100.root")             tmpfileName = "GaisserPDF_LowEnergy.root";
    else if (tmpfileName == "GaisserPDF_4000-10000-100000-1000-10000-0-1.5708-100.root")     tmpfileName = "GaisserPDF_MidEnergy.root";
    else if (tmpfileName == "GaisserPDF_100000-500000-1e007-10000-100000-0-1.5708-100.root") tmpfileName = "GaisserPDF_HighEnergy.root";
    
    std::ostringstream pdfFilePath;
    pdfFilePath << fInputDir << tmpfileName;
    std::string fileName = pdfFilePath.str();
    std::cout << "File path; " << fileName << std::endl;

    if(fSetRead){
      struct stat buffer;
      fSetRead = stat(fileName.c_str(), &buffer) == 0; // Check if file exists already
      if(!fSetRead) std::cout << "WARNING- "+fileName+" does not exist." << std::endl;
      else{
	std::cout << "Reading PDF from file "+fileName << std::endl;
	m_File = new TFile(fileName.c_str()); // Open file
	if(m_File->IsZombie() || m_File->TestBit(TFile::kRecovered)){ // Check that file is not corrupted
	  std::cout << "WARNING- "+fileName+" is corrupted or cannot be read." << std::endl;
	  fSetRead = false;
	}
      }
    }
    
    if(fSetRead){ // If the file exists then read it....
      std::cout << "Now going to read file...." << std::endl;
      fSetWrite = false; // Don't want to write as already exists
      double thetalow  = fThetamin;
      m_thetaHist = (TH1D*) m_File->Get("pdf_theta");
      for(int i=1; i<=fThetaBins; i++){	
	std::ostringstream pdfEnergyHist;
	pdfEnergyHist << "pdf_energy_"<<i;
	std::string pdfEnergyHiststr = pdfEnergyHist.str();
	
	TH1* pdf_hist = (TH1D*) m_File->Get( pdfEnergyHiststr.c_str() ); // Get the bin
	m_PDFmap->insert(std::make_pair(-thetalow,pdf_hist));         //---- -ve theta for quicker sorting
	thetalow = double(i)*(fThetamax)/double(fThetaBins); //---- increment the value of theta
      }
    }      // ------------------------------------------
    else { // File doesn't exist so want to make it..... 
           // ------------------------------------------
      std::cout << "Generating a new muon flux PDF" << std::endl;
      if(fSetWrite){
	std::cout << "Writing to PDF to file "+fileName << std::endl;
	m_File = new TFile(fileName.c_str(),"RECREATE");
      }
    
      double dnbins_theta = double(fThetaBins);
      m_thetaHist = new TH1D("pdf_theta", "pdf_theta", fThetaBins, fThetamin, fThetamax);
      for(int i=1; i<=fThetaBins; i++){
	double di = i==0 ? 0.1 : double(i);
	double theta = di*(fThetamax)/dnbins_theta;
	double int_i = muonSpec->Integral(fEmin, fEmax, fThetamin, theta, fEpsilon);
	m_thetaHist->SetBinContent(i, int_i/TotalMuonFlux);
      }
      if(fSetWrite) m_thetaHist->Write();
      std::cout << "theta PDF complete... now making the energy PDFs (this will take a little longer)... " << std::endl;
      
      //---- now compute the energy pdf
      double thetalow  = fThetamin;
      for(int i=1; i<=fThetaBins; i++){
      
	double di = double(i);
	double theta = di*(fThetamax)/fThetaBins;
      
	//---- compute the total integral
	double int_tot = muonSpec->Integral(fEmin, fEmax, thetalow, theta, fEpsilon);    
      
	//---- compute pdf for the low energy
	int nbins = fEBinsLow;
	TH1* pdf_lowenergy = new TH1D("pdf_lowenergy", "pdf_lowenergy", nbins, fEmin, fEmid);
	double dnbins = double(nbins);
	for(int j=1; j<=nbins; j++){
	  double dj = double(j);
	  double int_j = muonSpec->Integral(fEmin, fEmin + dj*(fEmid-fEmin)/dnbins, thetalow, theta, fEpsilon);
	  //	std::std::cout << j << "(" << m_emin << " --> " << m_emin + dj*m_emid/dnbins << ") = " << int_j/int_tot << std::std::endl;
	  //	std::std::cout << j << "(" << m_emin << " --> " << m_emin + dj*(m_emid-m_emin)/dnbins << ") = " << int_j/int_tot << std::std::endl;
	  pdf_lowenergy->SetBinContent(j, int_j/int_tot);	
	}
      
	//---- compute pdf for the high energy
	nbins = fEBinsHigh;
	dnbins=double(nbins);
	TH1* pdf_highenergy = new TH1D("pdf_highenergy", "pdf_highenergy", nbins, fEmid, fEmax);
	for(int j=1; j<=nbins; j++){
	  double dj = double(j);
	  double int_j = muonSpec->Integral(fEmin, fEmid + dj*(fEmax-fEmid)/dnbins, thetalow, theta, fEpsilon);
	  //      std::cout << j << "(" << m_emin << " --> " << m_emid + dj*(m_emax-m_emid)/dnbins << ") = " << int_j/int_tot << std::endl;
	  pdf_highenergy->SetBinContent(j, int_j/int_tot);	
	}

	//---- now combine the two energy hists
	std::vector<double> vxbins  = GetBinning(pdf_lowenergy->GetXaxis(),false);
	std::vector<double> vxbins2 = GetBinning(pdf_highenergy->GetXaxis());
	vxbins.insert(vxbins.end(), vxbins2.begin(), vxbins2.end());
      
	int ibin = 0;
	double* xbins = new double[vxbins.size()];
	for(std::vector<double>::const_iterator binit=vxbins.begin(); binit!=vxbins.end(); binit++, ibin++) xbins[ibin]=(*binit);
	TH1* pdf_energy = new TH1D("pdf_energy", "pdf_energy", vxbins.size()-1, xbins);
	int ibin2 = 1;
	for(ibin = 1; ibin<=pdf_lowenergy->GetNbinsX();   ibin++, ibin2++){
	  double content = pdf_lowenergy->GetBinContent(ibin);
	  if(ibin == 1) content = content - 0.00001;
	  pdf_energy->Fill(pdf_energy->GetBinCenter(ibin2), content);
	}
	for(ibin = 1; ibin<=pdf_highenergy->GetNbinsX();  ibin++, ibin2++){
	  pdf_energy->Fill(pdf_energy->GetBinCenter(ibin2), pdf_highenergy->GetBinContent(ibin));
	}
      
	//---- and remove any negative bins
	std::ostringstream Clonestr;
	Clonestr << "pdf_energy_"<<i;
	TH1* pdf_energy_noneg = (TH1D*) pdf_energy->Clone( Clonestr.str().c_str() );
	pdf_energy_noneg->Reset();
      
	double PDF = 0.0;
	double lastPD = 0.0;
	int nSkip=0;
	for(ibin = 1; ibin<=pdf_energy->GetNbinsX(); ibin++){
	  double newPD    = pdf_energy->GetBinContent(ibin);
	  double probDiff = newPD - lastPD;      
	  if(probDiff<0){
	    if(ibin!=pdf_energy->GetNbinsX()){
	      nSkip++;
	      continue;
	    }
	    else probDiff = 0;
	  }
	  else PDF += probDiff;
	  if(PDF > 1) PDF = 1;
	  for(int iskip=0; iskip <= nSkip; iskip++) pdf_energy_noneg->Fill(pdf_energy_noneg->GetBinCenter(ibin-iskip), PDF);
	  nSkip=0;	
	  lastPD = newPD;
	}
	
         
	//---- add this hist, increment thetalow and delete unwanted TH1s
	if(fSetWrite) pdf_energy_noneg->Write();
	m_PDFmap->insert(std::make_pair(-thetalow,pdf_energy_noneg)); //---- -ve theta for quicker sorting
    
	//---- increment the value of theta
	thetalow = theta; 
      
	//---- free up memory from unwanted hists
	delete pdf_lowenergy;
	delete pdf_highenergy;
	delete pdf_energy;
    
	std::cout << "\r===> " << 100.0*double(i)/double(fThetaBins) << "% complete... " << std::endl;
      } // ThetaBins
      std::cout << "finished the energy pdfs." << std::endl;
    }//---- if(!m_doRead)

    delete muonSpec;
    return;
  } // Make PDF

  //_____________________________________________________________________________
  double GaisserParam::GaisserFlux(double e, double theta){
  
    double ct  = cos(theta);    
    double di;
    if(fSetParam){
      //    double gamma=2.77;   // LVD spectrum: spectral index
      //    double A=1.84*0.14;  // normalisation
      double gamma = 2.7;
      double A = 0.14;
      double rc = 1.e-4;   // fraction of prompt muons  
      double c1  = sqrt(1.-(1.-pow(ct,2.0))/pow( (1.+32./6370.) ,2.0));  // Earth curvature    
      double deltae = 2.06e-3 * (1030. / c1 - 120.);         // muon energy loss in atmosphere    
      double em = e + deltae/2.;
      double e1 = e + deltae;
      double pdec = pow( (120. / (1030. / c1)), (1.04 / c1 / em));   // muon decay  
      di=A*pow(e1,-gamma)*(1./(1.+1.1*e1*c1/115.) + 0.054/(1.+1.1*e1*c1/850.) + rc)*pdec;
    }
    else{
      double gamma=2.7; // spectral index
      double A=0.14;    // normalisation
      double C = 3.64;
      double gamma2 = 1.29;
      double ct_star = sqrt( (pow(ct,2) + 0.102573*0.102573 - 0.068287*pow(ct,0.958633) 
			      + 0.0407253*pow(ct,0.817285) )/(1.0 + 0.102573*0.102573 - 0.068287 + 0.0407253) );
      double eMod = e*(1.0+C/(e*pow(ct_star,gamma2)));
      di=A*pow(eMod,-gamma)*(1./(1.+1.1*e*ct_star/115.) + 0.054/(1.+1.1*e*ct_star/850.));
    }
  
    return di;
  } // GaisserFlux

  //______________________________________________________________________________
  double GaisserParam::GaisserMuonFlux_Integrand(Double_t *x, Double_t*){  

    //---- calculate the flux
    double flux = 2.0*M_PI*sin(x[1])*GaisserFlux(x[0],x[1]);

    return flux;
  } // MuonFluxIntegrand

  //__________________________________
  std::vector<double> GaisserParam::GetBinning(const TAxis* axis, bool finalEdge){
    std::vector<double> vbins;
    for(int ibin=1; ibin<=axis->GetNbins()+1; ibin++){
      if(ibin<=axis->GetNbins()) vbins.push_back(axis->GetBinLowEdge(ibin));
      else if(finalEdge) vbins.push_back(axis->GetBinLowEdge(ibin) + axis->GetBinWidth(ibin));
    }
    return vbins;
  } // Get Binning 
  
  //__________________________________
  void GaisserParam::ResetMap(){  
    if(m_thetaHist) delete m_thetaHist;
    if(m_PDFmap){
      for(dhist_Map_it mapit=m_PDFmap->begin(); mapit!=m_PDFmap->end(); mapit++){
	if(mapit->second) delete mapit->second;
      }
      m_PDFmap->clear();  
      delete m_PDFmap;
      std::cout << "Reset PDFmap and thetaHist..." << std::endl;
    }
  } // ResetMap

}//end namespace evgen

namespace evgen{

  DEFINE_ART_MODULE(GaisserParam)

}//end namespace evgen

#endif
////////////////////////////////////////////////////////////////////////
