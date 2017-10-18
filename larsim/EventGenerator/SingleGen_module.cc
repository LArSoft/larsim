////////////////////////////////////////////////////////////////////////
/// \file  SingleGen_module.cc
/// \brief Generator for cosmic-rays
///
/// Module designed to produce a set list of particles for a MC event
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_SINGLEGEN
#define EVGEN_SINGLEGEN

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <memory>
#include <iterator>
#include <vector>


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
#include "nutools/RandomUtils/NuRandomService.h"

// nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"

#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace simb { class MCTruth; }

namespace evgen {

  /// module to produce single or multiple specified particles in the detector
  class SingleGen : public art::EDProducer {

  public:
    explicit SingleGen(fhicl::ParameterSet const& pset);
    virtual ~SingleGen();

    // This is called for each event.
    void produce(art::Event& evt);
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    void SampleOne(unsigned int   i, 
		   simb::MCTruth &mct);        
    void SampleMany(simb::MCTruth &mct);        
    void Sample(simb::MCTruth &mct);        
    void printVecs(std::vector<std::string> const& list);
    bool PadVector(std::vector<double> &vec);      
    double SelectFromHist(const TH1 *h);
    void SelectFromHist(const TH2 *h, double &x, double &y);
    
    static const int kUNIF = 0;    
    static const int kGAUS = 1;    
    static const int kHIST = 2;    // histogram for distribution

    int                 fMode;           ///< Particle Selection Mode 
                                         ///< 0--generate a list of all particles, 
                                         ///< 1--generate a single particle selected randomly from the list
    bool                fPadOutVectors;  ///< Select to pad out configuration vectors if they are not of 
					 ///< of the same length as PDG  
                                         ///< false: don't pad out - all values need to specified
                                         ///< true: pad out - default values assumed and printed out
    std::vector<int>    fPDG;            ///< PDG code of particles to generate    
    std::vector<double> fP0;             ///< Central momentum (GeV/c) to generate    
    std::vector<double> fSigmaP;         ///< Variation in momenta (GeV/c)    
    int                 fPDist;          ///< How to distribute momenta (gaus or uniform)    
    std::vector<double> fX0;             ///< Central x position (cm) in world coordinates 
    std::vector<double> fY0;             ///< Central y position (cm) in world coordinates
    std::vector<double> fZ0;             ///< Central z position (cm) in world coordinates
    std::vector<double> fT0;             ///< Central t position (s) in world coordinates
    std::vector<double> fSigmaX;         ///< Variation in x position (cm)    
    std::vector<double> fSigmaY;         ///< Variation in y position (cm)    
    std::vector<double> fSigmaZ;         ///< Variation in z position (cm)    
    std::vector<double> fSigmaT;         ///< Variation in t position (s)    
    int                 fPosDist;        ///< How to distribute xyz (gaus, or uniform)        
    int                 fTDist;          ///< How to distribute t  (gaus, or uniform)        
    bool                fSingleVertex;   ///< if true - all particles produced at the same location        
    std::vector<double> fTheta0XZ;       ///< Angle in XZ plane (degrees)    
    std::vector<double> fTheta0YZ;       ///< Angle in YZ plane (degrees)    
    std::vector<double> fSigmaThetaXZ;   ///< Variation in angle in XZ plane    
    std::vector<double> fSigmaThetaYZ;   ///< Variation in angle in YZ plane    
    int                 fAngleDist;      ///< How to distribute angles (gaus, uniform)
    std::string fHistFileName;               ///< Filename containing histogram of momenta
    std::vector<std::string> fPHist;     ///< name of histogram of momenta
//    std::vector<std::string> fThetaPhiHist; ///< name of histogram for theta/phi distribution
    std::vector<std::string> fThetaXzYzHist;   ///< name of histogram for thetaxz/thetayz distribution

    std::vector<TH1*> hPHist ;           /// actual TH1 for momentum distributions
//    std::vector<TH2*> hThetaPhiHist ;       /// actual TH1 for theta distributions - Theta on x axis
    std::vector<TH2*> hThetaXzYzHist ;         /// actual TH2 for angle distributions - Xz on x axis . 
    // FYI - thetaxz and thetayz are related to standard polar angles as follows:
    // thetaxz = atan2(math.sin(theta) * cos(phi), cos(theta))
    // thetayz = asin(sin(theta) * sin(phi));

    
  };
}

namespace evgen{

  //____________________________________________________________________________
  SingleGen::SingleGen(fhicl::ParameterSet const& pset)
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
  SingleGen::~SingleGen()
  {
  }

  //____________________________________________________________________________
  void SingleGen::reconfigure(fhicl::ParameterSet const& p)
  {
    // do not put seed in reconfigure because we don't want to reset 
    // the seed midstream
    fPadOutVectors = p.get< bool                >("PadOutVectors");
    fMode          = p.get< int                 >("ParticleSelectionMode");
    fPDG           = p.get< std::vector<int>    >("PDG");
    fPDist         = p.get< int                 >("PDist", kHIST);
    fTDist         = p.get< int                 >("TDist", kHIST);
    fX0            = p.get< std::vector<double> >("X0");
    fY0            = p.get< std::vector<double> >("Y0");
    fZ0            = p.get< std::vector<double> >("Z0");
    fT0            = p.get< std::vector<double> >("T0");
    fSigmaX        = p.get< std::vector<double> >("SigmaX");
    fSigmaY        = p.get< std::vector<double> >("SigmaY");
    fSigmaZ        = p.get< std::vector<double> >("SigmaZ");
    fSigmaT        = p.get< std::vector<double> >("SigmaT");
    fPosDist       = p.get< int                 >("PosDist");
    fHistFileName  = p.get< std::string         >("HistogramFile","");
    fPHist       = p.get< std::vector<std::string> >("PHist",std::vector<std::string>());
    fP0          = p.get< std::vector<double> >("P0");
    fSigmaP      = p.get< std::vector<double> >("SigmaP");
    fAngleDist     = p.get< int                 >("AngleDist");
//    fThetaPhiHist   = p.get< std::vector<std::string> >("ThetaPhiHist");
    fThetaXzYzHist     = p.get< std::vector<std::string> >("ThetaXzYzHist",std::vector<std::string>());
    fTheta0XZ    = p.get< std::vector<double> >("Theta0XZ");
    fTheta0YZ    = p.get< std::vector<double> >("Theta0YZ");
    fSigmaThetaXZ= p.get< std::vector<double> >("SigmaThetaXZ");
    fSigmaThetaYZ= p.get< std::vector<double> >("SigmaThetaYZ");

    std::vector<std::string> vlist(15);
    vlist[0]  = "PDG";
    vlist[1]  = "P0";
    vlist[2]  = "SigmaP";
    vlist[3]  = "X0";
    vlist[4]  = "Y0";
    vlist[5]  = "Z0";
    vlist[6]  = "SigmaX";
    vlist[7]  = "SigmaY";
    vlist[8]  = "SigmaZ";
    vlist[9]  = "Theta0XZ";
    vlist[10] = "Theta0YZ";
    vlist[11] = "SigmaThetaXZ";
    vlist[12] = "SigmaThetaYZ";
    vlist[13] = "T0";
    vlist[14] = "SigmaT";

//    vlist[15] = "PHist";
//    vlist[16] = "ThetaHist";
//    vlist[17] = "PhiHist";
    
    // begin tests for multiple particle error possibilities  
    std::string list;
    if( !this->PadVector(fP0          ) ){ list.append(vlist[1].append(", \n")); }
    if( !this->PadVector(fSigmaP      ) ){ list.append(vlist[2].append(", \n")); }
    if( !this->PadVector(fX0          ) ){ list.append(vlist[3].append(", \n")); }
    if( !this->PadVector(fY0          ) ){ list.append(vlist[4].append(", \n")); }
    if( !this->PadVector(fZ0          ) ){ list.append(vlist[5].append(", \n")); }
    if( !this->PadVector(fSigmaX      ) ){ list.append(vlist[6].append(", \n")); }
    if( !this->PadVector(fSigmaY      ) ){ list.append(vlist[7].append(", \n")); }
    if( !this->PadVector(fSigmaZ      ) ){ list.append(vlist[8].append(", \n")); }
    if( !this->PadVector(fTheta0XZ    ) ){ list.append(vlist[9].append(", \n")); }
    if( !this->PadVector(fTheta0YZ    ) ){ list.append(vlist[10].append(", \n")); }
    if( !this->PadVector(fSigmaThetaXZ) ){ list.append(vlist[11].append(", \n")); }
    if( !this->PadVector(fSigmaThetaYZ) ){ list.append(vlist[12].append("  \n")); }
    if( !this->PadVector(fT0          ) ){ list.append(vlist[13].append(", \n")); }
    if( !this->PadVector(fSigmaT      ) ){ list.append(vlist[14].append(", \n")); }

    

    if(list.size() > 0)
      throw cet::exception("SingleGen") << "The "<< list 
					<< "\n vector(s) defined in the fhicl files has/have "
					<< "a different size than the PDG vector "
					<< "\n and it has (they have) more than one value, "
					<< "\n disallowing sensible padding "
					<< " and/or you have set fPadOutVectors to false. \n";
    
    if(fPDG.size() > 1 && fPadOutVectors) this->printVecs(vlist);

    // If needed, get histograms for momentum and angle distributions
    TFile* fHistFile = nullptr;
    if (!fHistFileName.empty()) {
      fHistFile = new TFile(fHistFileName.c_str());
      if (!fHistFile->IsOpen()) {
        throw art::Exception(art::errors::NotFound)
          << "Can't open ROOT file from 'HistogramFile': \"" << fHistFileName << "\".";
      }
    }
    
    //
    // deal with momentum distribution
    //
    switch (fPDist) {
      case kHIST:
        if (!fHistFile) {
          throw art::Exception(art::errors::Configuration)
            << "Momentum distribution requested from histogram, but there is no 'HistogramFile' specified.";
        }
        if (fPHist.empty()) {
          throw art::Exception(art::errors::Configuration)
            << "Momentum distribution requested from histogram, but there is no 'PHist' specified.";
        }
        if (fPHist.size() != fPDG.size()) {
          throw art::Exception(art::errors::Configuration)
            << fPHist.size() << " momentum histograms to describe " << fPDG.size() << " particle types...";
        }
        hPHist.reserve(fPDG.size());
        for (auto const& histName: fPHist) {
          TH1* pHist = dynamic_cast<TH1*>(fHistFile->Get(histName.c_str()));
          if (!pHist) {
            throw art::Exception(art::errors::NotFound)
             << "Failed to read momentum histogram '" << histName << "' from '" << fHistFile->GetPath() << "\'";
          }
          pHist->SetDirectory(nullptr); // make it independent of the input file
          hPHist.emplace_back(pHist);
        } // for
        break;
      default:
        if (!fPHist.empty()) {
          throw art::Exception(art::errors::Configuration)
            << "Momentum distribution histograms specified, but there is no request to use them ('PDist').";
        }
        break;
    } // switch(fPDist)
    
    switch (fAngleDist) {
      case kHIST:
        if (!fHistFile) {
          throw art::Exception(art::errors::Configuration)
            << "Direction distribution requested from histogram, but there is no 'HistogramFile' specified.";
        }
        if (fThetaXzYzHist.empty()) {
          throw art::Exception(art::errors::Configuration)
            << "Direction distribution requested from histogram, but there is no 'ThetaXzYzHist' specified.";
        }
        if (fThetaXzYzHist.size() != fPDG.size()) {
          throw art::Exception(art::errors::Configuration)
            << fThetaXzYzHist.size() << " direction histograms to describe " << fPDG.size() << " particle types...";
        }
        hThetaXzYzHist.reserve(fPDG.size());
        for (auto const& histName: fThetaXzYzHist) {
          TH2* pHist = dynamic_cast<TH2*>(fHistFile->Get(histName.c_str()));
          if (!pHist) {
            throw art::Exception(art::errors::NotFound)
             << "Failed to read direction histogram '" << histName << "' from '" << fHistFile->GetPath() << "\'";
          }
          pHist->SetDirectory(nullptr); // make it independent of the input file
          hThetaXzYzHist.emplace_back(pHist);
        } // for
      default:
        if (!fThetaXzYzHist.empty()) {
          throw art::Exception(art::errors::Configuration)
            << "Direction distribution histograms specified, but there is no request to use them ('AngleDist').";
        }
        break;
    } // switch(fAngleDist)
    
    delete fHistFile;
    
#if 0
    if (fThetaXzYzHist.size() ==0){fAngleDist=0;}
    if (fPHist.size() ==0){fPDist=0;}
    if (fPDist==kHIST){
      hPHist.reserve(fPDG.size());
      for (unsigned int i(0); i < fPDG.size(); ++i){
        hPHist.emplace_back( (TH1*)fHistFile->Get( fPHist[i].c_str() ) );
      }
    }
    if (fAngleDist==kHIST){
//      hThetaPhiHist.reserve(fPDG.size());
      hThetaXzYzHist.reserve(fPDG.size());
      for (unsigned int i(0); i < fPDG.size(); ++i){
//        hThetaPhiHist.emplace_back( (TH2*)fHistFile->Get( fThetaPhiHist[i].c_str() ));
        hThetaXzYzHist.emplace_back( (TH2*)fHistFile->Get( fThetaXzYzHist[i].c_str() ));
      }
    }
#endif // 0
    return;
  }

  //____________________________________________________________________________
  bool SingleGen::PadVector(std::vector<double> &vec)
  {
    // check if the vec has the same size as fPDG
    if( vec.size() != fPDG.size() ){
      // if not padding out the vectors always cause an 
      // exception to be thrown if the vector in question
      // is not the same size as the fPDG vector
      // the exception is thrown in the reconfigure method
      // that calls this one
      if     (!fPadOutVectors) return false;
      else if( fPadOutVectors){
	// if padding of vectors is desired but the vector in
	// question has more than one entry it isn't clear
	// what the padded values should be so cause
	// an exception
	if(vec.size() != 1) return false;

	// pad it out
	vec.resize(fPDG.size(), vec[0]);

      }// end if padding out vectors
    }// end if the vector size is not the same as fPDG

    return true;
  }

  //____________________________________________________________________________
  void SingleGen::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));

    run.put(std::move(runcol));

    return;
  }

  //____________________________________________________________________________
  void SingleGen::produce(art::Event& evt)
  {

    ///unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    Sample(truth);

    LOG_DEBUG("SingleGen") << truth;

    truthcol->push_back(truth);

    evt.put(std::move(truthcol));

    return;
  }

  //____________________________________________________________________________
  // Draw the type, momentum and position of a single particle from the 
  // FCIHL description
  void SingleGen::SampleOne(unsigned int i, simb::MCTruth &mct){

    // get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat   flat(engine);
    CLHEP::RandGaussQ gauss(engine);

    // Choose momentum
    double p = 0.0;
    double m = 0.0;
    if (fPDist == kGAUS) {
      p = gauss.fire(fP0[i], fSigmaP[i]);
    }
    else if (fPDist == kHIST){
      p = SelectFromHist(hPHist[i]);
    }
    else{// if (fPDist == kUNIF) {
      p = fP0[i] + fSigmaP[i]*(2.0*flat.fire()-1.0);
    }
//    else {std::cout << "do not understand the value of PDist!";}

    static TDatabasePDG  pdgt;
    TParticlePDG* pdgp = pdgt.GetParticle(fPDG[i]);
    if (pdgp) m = pdgp->Mass();
    
    // Choose position
    TVector3 x;
    if (fPosDist == kGAUS) {
      x[0] = gauss.fire(fX0[i], fSigmaX[i]);;
      x[1] = gauss.fire(fY0[i], fSigmaY[i]);
      x[2] = gauss.fire(fZ0[i], fSigmaZ[i]);
    }
    else {
      x[0] = fX0[i] + fSigmaX[i]*(2.0*flat.fire()-1.0);
      x[1] = fY0[i] + fSigmaY[i]*(2.0*flat.fire()-1.0);
      x[2] = fZ0[i] + fSigmaZ[i]*(2.0*flat.fire()-1.0);
    }

    double t = 0.;
    if(fTDist==kGAUS){
      t = gauss.fire(fT0[i], fSigmaT[i]);
    }
    else{
      t = fT0[i] + fSigmaT[i]*(2.0*flat.fire()-1.0);
    }

    TLorentzVector pos(x[0], x[1], x[2], t);
    
    // Choose angles
    double thxz = 0;
    double thyz = 0;
    
    double thyzrads = 0; 
    double thyzradsplussigma = 0;
    double thyzradsminussigma = 0;

    if (fAngleDist == kGAUS) {
      thxz = gauss.fire(fTheta0XZ[i], fSigmaThetaXZ[i]);
      thyz = gauss.fire(fTheta0YZ[i], fSigmaThetaYZ[i]);
    }
    else if (fAngleDist == kHIST){ // Select thetaxz and thetayz from histogram
      double thetaxz = 0;
      double thetayz = 0;
      SelectFromHist(hThetaXzYzHist[i], thetaxz, thetayz);
      thxz = (180./M_PI)*thetaxz;
      thyz = (180./M_PI)*thetayz;
    }
    else {
      
      // Choose angles flat in phase space, which is flat in theta_xz 
      // and flat in sin(theta_yz).
   
      thxz = fTheta0XZ[i] + fSigmaThetaXZ[i]*(2.0*flat.fire()-1.0);
     
      thyzrads = std::asin(std::sin((M_PI/180.)*(fTheta0YZ[i]))); //Taking asin of sin gives value between -Pi/2 and Pi/2 regardless of user input
      thyzradsplussigma = TMath::Min((thyzrads + ((M_PI/180.)*fabs(fSigmaThetaYZ[i]))), M_PI/2.);
      thyzradsminussigma = TMath::Max((thyzrads - ((M_PI/180.)*fabs(fSigmaThetaYZ[i]))), -M_PI/2.);
         
      //uncomment line to print angular variation info
      //std::cout << "Central angle: " << (180./M_PI)*thyzrads << " Max angle: " << (180./M_PI)*thyzradsplussigma << " Min angle: " << (180./M_PI)*thyzradsminussigma << std::endl; 

      double sinthyzmin = std::sin(thyzradsminussigma);
      double sinthyzmax = std::sin(thyzradsplussigma);
      double sinthyz = sinthyzmin + flat.fire() * (sinthyzmax - sinthyzmin);
      thyz = (180. / M_PI) * std::asin(sinthyz);
    }
    
    double thxzrad=thxz*M_PI/180.0;	
    double thyzrad=thyz*M_PI/180.0;

    TLorentzVector pvec(p*std::cos(thyzrad)*std::sin(thxzrad),
			p*std::sin(thyzrad),
			p*std::cos(thxzrad)*std::cos(thyzrad),
			std::sqrt(p*p+m*m));
 
    // set track id to -i as these are all primary particles and have id <= 0
    int trackid = -1*(i+1);
    std::string primary("primary");

    simb::MCParticle part(trackid, fPDG[i], primary);
    part.AddTrajectoryPoint(pos, pvec);

    //std::cout << "Px: " <<  pvec.Px() << " Py: " << pvec.Py() << " Pz: " << pvec.Pz() << std::endl;
    //std::cout << "x: " <<  pos.X() << " y: " << pos.Y() << " z: " << pos.Z() << " time: " << pos.T() << std::endl;
    //std::cout << "YZ Angle: " << (thyzrad * (180./M_PI)) << " XZ Angle: " << (thxzrad * (180./M_PI)) << std::endl; 
     
    mct.Add(part);
  }

  //____________________________________________________________________________
  // Draw the type, momentum and position for all particles from the 
  // FCIHL description.  Start positions will all match but momenta and angles drawn from
  // distributions defined in the fhicls
  void SingleGen::SampleMany(simb::MCTruth &mct){

    // get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat   flat(engine);
    CLHEP::RandGaussQ gauss(engine);

    // Choose position
    TVector3 x;
    if (fPosDist == kGAUS) {
      x[0] = gauss.fire(fX0[0], fSigmaX[0]);;
      x[1] = gauss.fire(fY0[0], fSigmaY[0]);
      x[2] = gauss.fire(fZ0[0], fSigmaZ[0]);
    }
    else {
      x[0] = fX0[0] + fSigmaX[0]*(2.0*flat.fire()-1.0);
      x[1] = fY0[0] + fSigmaY[0]*(2.0*flat.fire()-1.0);
      x[2] = fZ0[0] + fSigmaZ[0]*(2.0*flat.fire()-1.0);
    }

    double t = 0.;
    if(fTDist==kGAUS){
      t = gauss.fire(fT0[0], fSigmaT[0]);
    }
    else{
      t = fT0[0] + fSigmaT[0]*(2.0*flat.fire()-1.0);
    }

    TLorentzVector pos(x[0], x[1], x[2], t);
    
    // loop through particles and select momenta and angles
    for (unsigned int i(0); i<fPDG.size(); ++i){
      // Choose momentum
      double p = 0.0;
      double m = 0.0;
      if (fPDist == kGAUS) {
        p = gauss.fire(fP0[i], fSigmaP[i]);
      }
      else if (fPDist == kHIST){
        p = SelectFromHist(hPHist[i]);
      }
      else {
        p = fP0[i] + fSigmaP[i]*(2.0*flat.fire()-1.0);
      }
  
      static TDatabasePDG  pdgt;
      TParticlePDG* pdgp = pdgt.GetParticle(fPDG[i]);
      if (pdgp) m = pdgp->Mass();
     
  
      // Choose angles
      double thxz = 0;
      double thyz = 0;
      
      double thyzrads = 0; 
      double thyzradsplussigma = 0;
      double thyzradsminussigma = 0;
  
      if (fAngleDist == kGAUS) {
        thxz = gauss.fire(fTheta0XZ[i], fSigmaThetaXZ[i]);
        thyz = gauss.fire(fTheta0YZ[i], fSigmaThetaYZ[i]);
      }
      else if (fAngleDist == kHIST){
        double thetaxz = 0;
        double thetayz = 0;
        SelectFromHist(hThetaXzYzHist[i], thetaxz, thetayz);
        thxz = (180./M_PI)*thetaxz;
        thyz = (180./M_PI)*thetayz;
      }
      else {
        
        // Choose angles flat in phase space, which is flat in theta_xz 
        // and flat in sin(theta_yz).
     
        thxz = fTheta0XZ[i] + fSigmaThetaXZ[i]*(2.0*flat.fire()-1.0);
       
        thyzrads = std::asin(std::sin((M_PI/180.)*(fTheta0YZ[i]))); //Taking asin of sin gives value between -Pi/2 and Pi/2 regardless of user input
        thyzradsplussigma = TMath::Min((thyzrads + ((M_PI/180.)*fabs(fSigmaThetaYZ[i]))), M_PI/2.);
        thyzradsminussigma = TMath::Max((thyzrads - ((M_PI/180.)*fabs(fSigmaThetaYZ[i]))), -M_PI/2.);
           
        //uncomment line to print angular variation info
        //std::cout << "Central angle: " << (180./M_PI)*thyzrads << " Max angle: " << (180./M_PI)*thyzradsplussigma << " Min angle: " << (180./M_PI)*thyzradsminussigma << std::endl; 
  
        double sinthyzmin = std::sin(thyzradsminussigma);
        double sinthyzmax = std::sin(thyzradsplussigma);
        double sinthyz = sinthyzmin + flat.fire() * (sinthyzmax - sinthyzmin);
        thyz = (180. / M_PI) * std::asin(sinthyz);
      }
      
      double thxzrad=thxz*M_PI/180.0;	
      double thyzrad=thyz*M_PI/180.0;
  
      TLorentzVector pvec(p*std::cos(thyzrad)*std::sin(thxzrad),
  			p*std::sin(thyzrad),
  			p*std::cos(thxzrad)*std::cos(thyzrad),
  			std::sqrt(p*p+m*m));
   
      // set track id to -i as these are all primary particles and have id <= 0
      int trackid = -1*(i+1);
      std::string primary("primary");
  
      simb::MCParticle part(trackid, fPDG[i], primary);
      part.AddTrajectoryPoint(pos, pvec);
  
      //std::cout << "Px: " <<  pvec.Px() << " Py: " << pvec.Py() << " Pz: " << pvec.Pz() << std::endl;
      //std::cout << "x: " <<  pos.X() << " y: " << pos.Y() << " z: " << pos.Z() << " time: " << pos.T() << std::endl;
      //std::cout << "YZ Angle: " << (thyzrad * (180./M_PI)) << " XZ Angle: " << (thxzrad * (180./M_PI)) << std::endl; 
      mct.Add(part);
    }
  }


  //____________________________________________________________________________
  void SingleGen::Sample(simb::MCTruth &mct) 
  {

    switch (fMode) {
    case 0: // List generation mode: every event will have one of each
	    // particle species in the fPDG array
        if (fSingleVertex){
          SampleMany(mct);
        }
        else{
          for (unsigned int i=0; i<fPDG.size(); ++i) {
            SampleOne(i,mct);
          }//end loop over particles
        }
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
      mf::LogWarning("UnrecognizeOption") << "SingleGen does not recognize ParticleSelectionMode "
					  << fMode;
      break;
    } // switch on fMode

    return;
  }

  //____________________________________________________________________________
  void SingleGen::printVecs(std::vector<std::string> const& list)
  {
 
    mf::LogInfo("SingleGen") << " You are using vector values for SingleGen configuration.\n   " 
			     << " Some of the configuration vectors may have been padded out ,"
			     << " because they (weren't) as long as the pdg vector"
			     << " in your configuration. \n"
			     << " The new input particle configuration is:\n" ;

    std::string values;
    for(size_t i = 0; i <=1; ++i){// list.size(); ++i){

      values.append(list[i]);
      values.append(": [ ");      
      
      for(size_t e = 0; e < fPDG.size(); ++e){
        std::stringstream buf;
        buf.width(10);
	if(i == 0 ) buf << fPDG[e]          << ", ";
	buf.precision(5);
	if(i == 1 ) buf << fP0[e]           << ", ";
	if(i == 2 ) buf << fSigmaP[e] 	    << ", ";
	if(i == 3 ) buf << fX0[e]     	    << ", ";
	if(i == 4 ) buf << fY0[e]     	    << ", ";
	if(i == 5 ) buf << fZ0[e]	    << ", ";
	if(i == 6 ) buf << fSigmaX[e] 	    << ", ";
	if(i == 7 ) buf << fSigmaY[e] 	    << ", ";
	if(i == 8 ) buf << fSigmaZ[e] 	    << ", ";
	if(i == 9 ) buf << fTheta0XZ[e]     << ", ";
	if(i == 10) buf << fTheta0YZ[e]     << ", ";
	if(i == 11) buf << fSigmaThetaXZ[e] << ", ";
	if(i == 12) buf << fSigmaThetaYZ[e] << ", ";
	if(i == 13) buf << fT0[e]           << ", ";
	if(i == 14) buf << fSigmaT[e]       << ", ";
        values.append(buf.str());
      }

      values.erase(values.find_last_of(","));
      values.append(" ] \n");

    }// end loop over vector names in list

    mf::LogInfo("SingleGen") << values;

    return;
  }
  
  
  //____________________________________________________________________________
  double SingleGen::SelectFromHist(const TH1 *h) // select from a 1D histogram
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat   flat(engine);
    
    double throw_value = h->Integral() * flat.fire();
    double cum_value(0);
    for (int i(0); i < h->GetNbinsX()+1; ++i){
      cum_value += h->GetBinContent(i);
      if (throw_value < cum_value){
        return flat.fire()*h->GetBinWidth(i) + h->GetBinLowEdge(i);
      }
    }
    return throw_value; // for some reason we've gone through all bins and failed?
  }
  //____________________________________________________________________________
  void SingleGen::SelectFromHist(const TH2 *h, double &x, double &y) // select from a 2D histogram
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat   flat(engine);
    
    double throw_value = h->Integral() * flat.fire();
    double cum_value(0);
    for (int i(0); i < h->GetNbinsX()+1; ++i){
      for (int j(0); j < h->GetNbinsY()+1; ++j){
        cum_value += h->GetBinContent(i, j);
        if (throw_value < cum_value){
          x = flat.fire()*h->GetXaxis()->GetBinWidth(i) + h->GetXaxis()->GetBinLowEdge(i);
          y = flat.fire()*h->GetYaxis()->GetBinWidth(j) + h->GetYaxis()->GetBinLowEdge(j);
          return;
        }
      }
    }
    return; // for some reason we've gone through all bins and failed?
  }
  //____________________________________________________________________________


}//end namespace evgen

namespace evgen{

  DEFINE_ART_MODULE(SingleGen)

}//end namespace evgen

#endif
////////////////////////////////////////////////////////////////////////
