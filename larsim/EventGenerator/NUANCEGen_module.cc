////////////////////////////////////////////////////////////////////////
// $Id: NUANCEGen_module.cc,v 1.4 2010/04/27 19:48:46 brebel Exp $
//
//
// NUANCE neutrino event generator
//
// brebel@fnal.gov
// saima@ksu.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_NUANCEGEN_H
#define EVGEN_NUANCEGEN_H

#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <unistd.h>
#include <stdio.h>
#include <fstream>

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TStopwatch.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"

// LArSoft includes
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCNeutrino.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/SummaryData/RunData.h"

// class TH1F;
// class TH2F;

// namespace simb { class MCTruth;     }

namespace evgen {
  /// A module to check the results from the Monte Carlo generator
  class NUANCEGen : public art::EDProducer {
  public:
    explicit NUANCEGen(fhicl::ParameterSet const &pset);
    virtual ~NUANCEGen();                        

    void produce(art::Event& evt);  
    void beginJob();
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);
    void endJob();

  private:

        std::string ParticleStatus(int StatusCode);
        std::string ReactionChannel(int ccnc,int mode);
    
        void FillHistograms(simb::MCTruth mc);

        std::string         fNuanceFile;
        double              fBeamVerticalAngle;
	std::ifstream      *fEventFile;
	TStopwatch          fStopwatch;      ///keep track of how long it takes to run the job
	
	
	std::string fNUANCEModuleLabel;
	
	TH1F* fGenerated[6];  ///< Spectra as generated
	
	TH1F* fVertexX;    ///< vertex location of generated events in x
	TH1F* fVertexY;    ///< vertex location of generated events in y
	TH1F* fVertexZ;    ///< vertex location of generated events in z
	
	TH2F* fVertexXY;   ///< vertex location in xy
	TH2F* fVertexXZ;   ///< vertex location in xz
	TH2F* fVertexYZ;   ///< vertex location in yz
	
	TH1F* fDCosX;      ///< direction cosine in x
	TH1F* fDCosY;      ///< direction cosine in y
	TH1F* fDCosZ;      ///< direction cosine in z
	
	TH1F* fMuMomentum; ///< momentum of outgoing muons
	TH1F* fMuDCosX;    ///< direction cosine of outgoing mu in x
	TH1F* fMuDCosY;    ///< direction cosine of outgoing mu in y
	TH1F* fMuDCosZ;    ///< direction cosine of outgoing mu in z
	
	TH1F* fEMomentum;  ///< momentum of outgoing electrons
	TH1F* fEDCosX;     ///< direction cosine of outgoing e in x
	TH1F* fEDCosY;     ///< direction cosine of outgoing e in y
	TH1F* fEDCosZ;     ///< direction cosine of outgoing e in z
	
	TH1F* fCCMode;      ///< CC interaction mode
	TH1F* fNCMode;      ///< CC interaction mode
	
	TH1F* fDeltaE;     ///< difference in neutrino energy from MCTruth::Enu() vs TParticle
	TH1F* fECons;      ///< histogram to determine if energy is conserved in the event
	
  };
}

namespace evgen{

  //____________________________________________________________________________
  NUANCEGen::NUANCEGen(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset); 
    fStopwatch.Start();

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();

    fEventFile = new ifstream(fNuanceFile.c_str());
   }

  //____________________________________________________________________________
  NUANCEGen::~NUANCEGen()
  {  
    fStopwatch.Stop();
  }

  //____________________________________________________________________________
  void NUANCEGen::reconfigure(fhicl::ParameterSet const& p)
  {
    fNuanceFile          =(p.get< std::string         >("NuanceFile"));
    fBeamVerticalAngle   =(p.get< double              >("BeamVerticalAngle"));
    return;
  }
//___________________________________________________________________________
  
  void NUANCEGen::beginJob(){
   
    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    fGenerated[0] = tfs->make<TH1F>("fGenerated_necc","",  100, 0.0, 20.0);
    fGenerated[1] = tfs->make<TH1F>("fGenerated_nebcc","", 100, 0.0, 20.0);
    fGenerated[2] = tfs->make<TH1F>("fGenerated_nmcc","",  100, 0.0, 20.0);
    fGenerated[3] = tfs->make<TH1F>("fGenerated_nmbcc","", 100, 0.0, 20.0);
    fGenerated[4] = tfs->make<TH1F>("fGenerated_nnc","",   100, 0.0, 20.0);
    fGenerated[5] = tfs->make<TH1F>("fGenerated_nbnc","",  100, 0.0, 20.0);
    
    fDCosX = tfs->make<TH1F>("fDCosX", ";dx/ds", 200, -1., 1.);
    fDCosY = tfs->make<TH1F>("fDCosY", ";dy/ds", 200, -1., 1.);
    fDCosZ = tfs->make<TH1F>("fDCosZ", ";dz/ds", 200, -1., 1.);

    fMuMomentum = tfs->make<TH1F>("fMuMomentum", ";p_{#mu} (GeV/c)", 500, 0., 50.);
    fMuDCosX    = tfs->make<TH1F>("fMuDCosX", ";dx/ds;", 200, -1., 1.);
    fMuDCosY    = tfs->make<TH1F>("fMuDCosY", ";dy/ds;", 200, -1., 1.);
    fMuDCosZ    = tfs->make<TH1F>("fMuDCosZ", ";dz/ds;", 200, -1., 1.);

    fEMomentum  = tfs->make<TH1F>("fEMomentum", ";p_{e} (GeV/c)", 500, 0., 50.);
    fEDCosX     = tfs->make<TH1F>("fEDCosX", ";dx/ds;", 200, -1., 1.);
    fEDCosY     = tfs->make<TH1F>("fEDCosY", ";dy/ds;", 200, -1., 1.);
    fEDCosZ     = tfs->make<TH1F>("fEDCosZ", ";dz/ds;", 200, -1., 1.);

    fCCMode = tfs->make<TH1F>("fCCMode", ";CC Interaction Mode;", 5, 0., 5.);
    fCCMode->GetXaxis()->SetBinLabel(1, "QE");
    fCCMode->GetXaxis()->SetBinLabel(2, "Res");
    fCCMode->GetXaxis()->SetBinLabel(3, "DIS");
    fCCMode->GetXaxis()->SetBinLabel(4, "Coh");
    fCCMode->GetXaxis()->SetBinLabel(5, "kInverseMuDecay");
    fCCMode->GetXaxis()->CenterLabels();

    fNCMode = tfs->make<TH1F>("fNCMode", ";NC Interaction Mode;", 5, 0., 5.);
    fNCMode->GetXaxis()->SetBinLabel(1, "QE");
    fNCMode->GetXaxis()->SetBinLabel(2, "Res");
    fNCMode->GetXaxis()->SetBinLabel(3, "DIS");
    fNCMode->GetXaxis()->SetBinLabel(4, "Coh");
    fNCMode->GetXaxis()->SetBinLabel(5, "kNuElectronElastic");
    fNCMode->GetXaxis()->CenterLabels();

    //fDeltaE = tfs->make<TH1F>("fDeltaE", ";#Delta E_{#nu} (GeV);", 200, -1., 1.); 
    fECons  = tfs->make<TH1F>("fECons", ";#Delta E(#nu,lepton);", 500, -5., 5.);

    art::ServiceHandle<geo::Geometry> geo;
    double x = 2.1*geo->DetHalfWidth();
    double y = 2.1*geo->DetHalfHeight();
    double z = 2.*geo->DetLength();
    int xdiv = TMath::Nint(2*x/5.);
    int ydiv = TMath::Nint(2*y/5.);
    int zdiv = TMath::Nint(2*z/5.);
    
    fVertexX = tfs->make<TH1F>("fVertexX", ";x (cm)", xdiv,  -x, x);
    fVertexY = tfs->make<TH1F>("fVertexY", ";y (cm)", ydiv,  -y, y);
    fVertexZ = tfs->make<TH1F>("fVertexZ", ";z (cm)", zdiv, -0.2*z, z);
    
    fVertexXY = tfs->make<TH2F>("fVertexXY", ";x (cm);y (cm)", xdiv,     -x, x, ydiv, -y, y);
    fVertexXZ = tfs->make<TH2F>("fVertexXZ", ";z (cm);x (cm)", zdiv, -0.2*z, z, xdiv, -x, x);
    fVertexYZ = tfs->make<TH2F>("fVertexYZ", ";z (cm);y (cm)", zdiv, -0.2*z, z, ydiv, -y, y);

  }

  //____________________________________________________________________________
  void NUANCEGen::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));

    run.put(std::move(runcol));

    return;
  }

  //____________________________________________________________________________
  void NUANCEGen::endJob()
  {
    fEventFile->close();
  }
  //____________________________________________________________________________
  void NUANCEGen::produce(art::Event& evt)
  {

    std::cout << std::endl;
    std::cout<<"------------------------------------------------------------------------------"<<std::endl;
//  std::cout << "run    : " << evt.Header().Run() << std::endl;
//  std::cout << "subrun : " << evt.Header().Subrun() << std::endl;
//  std::cout << "event  : " << evt.Header().Event() << std::endl;
    std::cout << "event  : " << evt.id().event() << std::endl;  
    
    double X = 0.; // vertex position from Nuance
    double Y = 0.; // vertex position from Nuance
    double Z = 0.; // vertex position from Nuance
    double PDGCODE = -9999.; 
    double CHANNEL = -9999.;
    int channel = -9999;
    double energy = 0.; // in MeV from Nuance
    double cosx = 0.;
    double cosy = 0.;
    double cosz = 0.;
    std::string name, k, dollar;
    int partnumber = 0;
    
    int trackid = -1; // set track id to -i as these are all primary particles and have id <= 0
    std::string primary("primary");
    int FirstMother = -1;
    double Mass = -9999;
    int Status = -9999;
    
    
    int ccnc = -9999;
    int mode = -9999;
    int targetnucleusPdg = -9999;
    int hitquarkPdg = -9999;

    double P; // momentum of MCParticle in GeV/c
    TLorentzVector Neutrino;
    TLorentzVector Lepton;
    TLorentzVector Target;
    TLorentzVector q;
    TLorentzVector Hadron4mom;

    double InvariantMass = -9999;
    double x = -9999;
    double y = -9999; 
    double Q2 = -9999;

    int Tpdg = 0;  // for target 
    double Tmass = 0;
    int Tstatus = 11;
    double Tcosx = 0., Tcosy = 0., Tcosz = 0., Tenergy = 0.;
    TLorentzVector Tpos;
    double M = 0;

    // rotated direction cosines
    double CosX = 0;
    double CosY = 0;
    double CosZ = 0;
    
        
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
    simb::MCTruth truth;
    
    if(!fEventFile->good())
      std::cout << "NuanceFile: Problem reading nuance file" << std::endl; 
    
    while(getline(*fEventFile,k)){
      std::istringstream in;
      in.clear();
      in.str(k);
      
      in>>dollar>>name>>PDGCODE>>energy>>cosx>>cosy>>cosz>>partnumber;
      //std::cout<<std::setprecision(9)<<dollar<<"  "<<name<<"  "<<PDGCODE<<"  "<<energy<<"  "<<cosx<<" "<<cosy<<"  "<<cosz<<"  "<<partnumber<<std::endl;
      
      //get the nuance channel number
      if(name == "nuance"){
	CHANNEL = PDGCODE;
	channel = int (CHANNEL);
	channel = abs(channel)+ simb::kNuanceOffset;
	//std::cout << "channel = " << channel << std::endl;
	
	//set the interaction type; CC or NC       
	if ( abs(channel) == simb::kCCQE 
	     || ( abs(channel) >= simb::kResCCNuProtonPiPlus && abs(channel) <= simb::kResCCNuNeutronPiPlus ) 
	     || ( abs(channel) >= simb::kResCCNuBarNeutronPiMinus && abs(channel) <= simb::kResCCNuBarProtonPiMinus ) 
	     || ( abs(channel) >=  simb::kResCCNuDeltaPlusPiPlus && abs(channel) <= simb::kResCCNuDelta2PlusPiMinus ) 
	     || ( abs(channel) >= simb::kResCCNuBarDelta0PiMinus && abs(channel) <= simb::kResCCNuBarDeltaMinusPiPlus ) 
	     || ( abs(channel) >= simb::kResCCNuProtonRhoPlus && abs(channel) <= simb::kResCCNuNeutronRhoPlus ) 
	     || ( abs(channel) >= simb::kResCCNuBarNeutronRhoMinus && abs(channel) <= simb::kResCCNuBarNeutronRho0 ) 
	     || ( abs(channel) >= simb::kResCCNuSigmaPlusKaonPlus && abs(channel) <= simb::kResCCNuSigmaPlusKaon0 )
	     || ( abs(channel) >= simb::kResCCNuBarSigmaMinusKaon0 && abs(channel) <= simb::kResCCNuBarSigma0Kaon0 ) 
	     || abs(channel) == simb::kResCCNuProtonEta 
	     || abs(channel) == simb::kResCCNuBarNeutronEta 
	     || abs(channel) == simb::kResCCNuKaonPlusLambda0
	     || abs(channel) == simb::kResCCNuBarKaon0Lambda0
	     || ( abs(channel) >= simb::kResCCNuProtonPiPlusPiMinus && abs(channel) <=  simb::kResCCNuProtonPi0Pi0 ) 
	     || ( abs(channel) >= simb::kResCCNuBarNeutronPiPlusPiMinus && abs(channel) <= simb::kResCCNuBarNeutronPi0Pi0 ) 
	     || abs(channel) == simb::kCCDIS || abs(channel) == simb::kCCQEHyperon 
	     || abs(channel) == simb::kCCCOH || abs(channel) == simb::kInverseMuDecay )
	  ccnc = simb::kCC;
	  	
	else if ( abs(channel) !=  simb::kUnUsed1 && abs(channel) !=  simb::kUnUsed2 )
	  ccnc = simb::kNC;
	  
	//set the interaction mode; QE, Res, DIS, Coh, kNuElectronElastic, kInverseMuDecay 
	if ( abs(channel) == simb::kCCQE || abs(channel) ==  simb::kNCQE || abs(channel) == simb::kCCQEHyperon )
	  mode = simb::kQE; 
	else if ( abs(channel) >= simb::kResCCNuProtonPiPlus && abs(channel) <= simb::kResCCNuBarProtonPi0Pi0 )
	  mode = simb::kRes;
	else if ( abs(channel) ==  simb::kCCDIS || abs(channel) ==  simb::kNCDIS )
	  mode = simb::kDIS;
	else if ( abs(channel) == simb::kNCCOH || abs(channel) == simb::kCCCOH )
	  mode = simb::kCoh;
	else if ( abs(channel) == simb::kNuElectronElastic )
	  mode = 4;
	else if ( abs(channel) == simb::kInverseMuDecay )
	  mode = 5;

      } // if(name == "nuance")


      //get the nuance vertex position:
      if(name == "vertex"){
	X = PDGCODE;
	Y = energy;
	Z = cosx;
	//std::cout << "vertex from nuance = " << X << "  " << Y << "  " << Z << std::endl;
      }

	double PI = 3.14159265;
	
	CosX = cosx;
	CosY = (TMath::Cos(fBeamVerticalAngle*PI/180)*cosy)+(TMath::Sin(fBeamVerticalAngle*PI/180)*cosz);
	CosZ = (-TMath::Sin(fBeamVerticalAngle*PI/180)*cosy)+(TMath::Cos(fBeamVerticalAngle*PI/180)*cosz);

	//get the target info
	if(name == "track" && (PDGCODE == 2212 || PDGCODE == 2112 || PDGCODE == 18040 || PDGCODE == 11) && partnumber == -1){
	  Tpdg = int (PDGCODE);
	
	  if ( PDGCODE == 18040 ){
	    Tmass = 37.5593438; // GeV
	  }	

	else {	
	  const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
	  const TParticlePDG* definition = databasePDG->GetParticle(Tpdg);
	  Tmass = definition->Mass();
	}
	  
	  Tenergy = energy;
	  Tcosx = CosX;
	  Tcosy = CosY;
	  Tcosz = CosZ;
	}
	
      //get mcparticles other than target
	if(name == "track" && ( 
			       partnumber == 0 || 
			       (partnumber == -1 && (PDGCODE != 2212 && PDGCODE != 2112 && PDGCODE != 18040 && PDGCODE != 11))
			       )
	   ){ 
	  
	  int pdgcode = int (PDGCODE); 
	  //std::cout << "pdgcode = " << pdgcode << std::endl;	
	  
	  if ( PDGCODE == 18040 )
	    Mass = 37.5593438; // GeV
	
	  else {
	    const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
	    const TParticlePDG* definition = databasePDG->GetParticle(pdgcode);
	    Mass = definition->Mass(); // GeV
	  }
	  
	  if(partnumber == -1)
	    Status = 0;
	  
	  if(partnumber == 0)
	    Status = 1;
	  
	  simb::MCParticle mcpart(trackid,
				  pdgcode,
				  primary,		    
				  FirstMother,
				  Mass,
				  Status
				  );
	  
	  P = std::sqrt(pow(energy/1000,2.) - pow(Mass,2.)); // GeV/c
	  //std::cout << "Momentum = " << P << std::endl;
	  
	  art::ServiceHandle<geo::Geometry> geo;
	  
	  double X0 = X + geo->DetHalfWidth();
	  double Y0 = Y;
	  double Z0 = Z + 0.5*geo->DetLength();
	  
	  TLorentzVector pos(X0, Y0, Z0, 0);
	  Tpos = pos; // for target
	  
	  TLorentzVector mom(CosX*P, CosY*P, CosZ*P, energy/1000);
	  
	  mcpart.AddTrajectoryPoint(pos,mom);
	  truth.Add(mcpart);
	  
	
	if(name == "track" && (abs(pdgcode) == 14 || abs(pdgcode) == 12) && partnumber == -1)
	  Neutrino.SetPxPyPzE(CosX*P, CosY*P, CosZ*P, energy/1000);
	
	if(name == "track" && (abs(pdgcode) == 13 || abs(pdgcode) == 11 || abs(pdgcode) == 14 || abs(pdgcode) == 12) && partnumber == 0)
	  Lepton.SetPxPyPzE(CosX*P, CosY*P, CosZ*P, energy/1000);

	}// loop over particles in an event. (excluding target)
      
	if(name == "end"){
	  break;  
	}
	
    } // end while loop

	    
    /////////////////////////////////
    
    //adding the target to mcparticle
    simb::MCParticle mcpart(trackid,
			    Tpdg,
			    primary,		    
			    FirstMother,
			    Tmass,
			    Tstatus
			    );
        
    TLorentzVector Tmom;
    Tmom.SetPxPyPzE(Tcosx, Tcosy, Tcosz, Tenergy/1000); // target momentum coordinates make a unit vector, may be modified later, does not effect the event simulation 

    mcpart.AddTrajectoryPoint(Tpos,Tmom); 
    truth.Add(mcpart);
    
//     std::cout<<"Neutrino 4-Momentum = "
// 	     <<Neutrino.Px()<<" "
// 	     << Neutrino.Py()<<" "
// 	     <<Neutrino.Pz()<<" "
// 	     <<Neutrino.E()<<std::endl;
    
//     std::cout << "Target 4-Momentum = " 
// 	      << Target.Px() << " " 
// 	      << Target.Py() << " " 
// 	      << Target.Pz()<< " " 
// 	      << Target.E() << std::endl;

//     std::cout << "4-momentum of Lepton = " 
// 	      << Lepton.Px() << " " 
// 	      << Lepton.Py() << " " 
// 	      << Lepton.Pz() << " " 
// 	      << Lepton.E() << std::endl; 
    
    M = Tmass;
    q = Neutrino - Lepton;
    Q2 = -(q*q);

    double v = q.E();
    x = Q2/ (2*M*v);
    y = v/ Neutrino.E();
    double W2 = M*M + 2*M*v - Q2;
    InvariantMass = TMath::Sqrt(TMath::Max(0.,W2));

    if(mode == simb::kCoh){
      x = -1;
      y = -1;
      InvariantMass = -1;
    }
    
    truth.SetOrigin(simb::kBeamNeutrino);
    truth.SetNeutrino(ccnc, mode, channel,
		      targetnucleusPdg, 
		      Tpdg, 
		      hitquarkPdg,
		      //InvariantMass, x, y, Q2
		      InvariantMass, x, y, Q2
		      );
    
    std::cout << truth.GetNeutrino() << std::endl;
     
     truthcol->push_back(truth);
     FillHistograms(truth);  
     evt.put(std::move(truthcol));
    
    return;
  }
  
//   //......................................................................
  std::string NUANCEGen::ParticleStatus(int StatusCode)
  {
    int code = StatusCode;
    std::string ParticleStatusName;

    switch(code)
      {
      case 0:
	ParticleStatusName = "kIStInitialState";
	break;
      case 1:
	ParticleStatusName = "kIStFinalState";
	break;
      case 11:
	ParticleStatusName = "kIStNucleonTarget";
	break;
      default:
	ParticleStatusName = "Status Unknown";
      }
    return ParticleStatusName;
  }


//   //......................................................................
  std::string NUANCEGen::ReactionChannel(int ccnc,int mode)
  {
    std::string ReactionChannelName=" ";

    if(ccnc==0)
      ReactionChannelName = "kCC";
    else if(ccnc==1)
      ReactionChannelName = "kNC";
    else std::cout<<"Current mode unknown!! "<<std::endl;

    if(mode==0)
      ReactionChannelName += "_kQE";
    else if(mode==1)
      ReactionChannelName += "_kRes";
    else if(mode==2)
      ReactionChannelName += "_kDIS";
    else if(mode==3)
      ReactionChannelName += "_kCoh";
    else if(mode==4)
      ReactionChannelName += "_kNuElectronElastic";
    else if(mode==5)
      ReactionChannelName += "_kInverseMuDecay";
    else std::cout<<"interaction mode unknown!! "<<std::endl;

    return ReactionChannelName;
  }

//   //......................................................................
  void NUANCEGen::FillHistograms(simb::MCTruth mc)
  {
    // Decide which histograms to put the spectrum in
    int id = -1;
    if (mc.GetNeutrino().CCNC()==simb::kCC) {
      fCCMode->Fill(mc.GetNeutrino().Mode());
      if      (mc.GetNeutrino().Nu().PdgCode() ==  12) id = 0;
      else if (mc.GetNeutrino().Nu().PdgCode() == -12) id = 1;
      else if (mc.GetNeutrino().Nu().PdgCode() ==  14) id = 2;
      else if (mc.GetNeutrino().Nu().PdgCode() == -14) id = 3;
      else return;
    }
    else {
      fNCMode->Fill(mc.GetNeutrino().Mode());
      if (mc.GetNeutrino().Nu().PdgCode() > 0) id = 4;
      else                                     id = 5;
    }
    if (id==-1) abort();
    
    // Fill the specta histograms
    fGenerated[id]->Fill(mc.GetNeutrino().Nu().E() );
      
    //< fill the vertex histograms from the neutrino - that is always 
    //< particle 0 in the list
    simb::MCNeutrino       mcnu = mc.GetNeutrino();
    const simb::MCParticle nu   = mcnu.Nu();
    
    fVertexX->Fill(nu.Vx());
    fVertexY->Fill(nu.Vy());
    fVertexZ->Fill(nu.Vz());
    
    fVertexXY->Fill(nu.Vx(), nu.Vy());
    fVertexXZ->Fill(nu.Vz(), nu.Vx());
    fVertexYZ->Fill(nu.Vz(), nu.Vy());
    
    double mom = nu.P();
    if(std::abs(mom) > 0.){
      fDCosX->Fill(nu.Px()/mom);
      fDCosY->Fill(nu.Py()/mom);
      fDCosZ->Fill(nu.Pz()/mom);
    }


//     LOG_DEBUG("GENIEInteractionInformation") 
//       << std::endl
//       << "REACTION:  " << ReactionChannel(mc.GetNeutrino().CCNC(),mc.GetNeutrino().Mode()) 
//       << std::endl
//       << "-----------> Particles in the Stack = " << mc.NParticles() << std::endl
//       << std::setiosflags(std::ios::left) 
//       << std::setw(20) << "PARTICLE"
//       << std::setiosflags(std::ios::left) 
//       << std::setw(32) << "STATUS"
//       << std::setw(18) << "E (GeV)"
//       << std::setw(18) << "m (GeV/c2)"
//       << std::setw(18) << "Ek (GeV)"
//       << std::endl << std::endl;

//     const TDatabasePDG* databasePDG = TDatabasePDG::Instance();

//     // Loop over the particle stack for this event 
//     for(int i = 0; i < mc.NParticles(); ++i){
//       simb::MCParticle part(mc.GetParticle(i));
//       std::string name = databasePDG->GetParticle(part.PdgCode())->GetName();
//       int code = part.StatusCode();
//       std::string status = ParticleStatus(code);
//       double mass = part.Mass();
//       double energy = part.E(); 
//       double Ek = (energy-mass); // Kinetic Energy (GeV)
//       if(status=="kIStFinalStB4Interactions")
// 	LOG_DEBUG("GENIEFinalState")
// 	  << std::setiosflags(std::ios::left) << std::setw(20) << name
// 	  << std::setiosflags(std::ios::left) << std::setw(32) <<status
// 	  << std::setw(18)<< energy
// 	  << std::setw(18)<< mass
// 	  << std::setw(18)<< Ek <<std::endl;
//       else 
// 	LOG_DEBUG("GENIEFinalState") 
// 	  << std::setiosflags(std::ios::left) << std::setw(20) << name
// 	  << std::setiosflags(std::ios::left) << std::setw(32) << status
// 	  << std::setw(18) << energy
// 	  << std::setw(18) << mass <<std::endl; 

    std::cout << "REACTION:  " << ReactionChannel(mc.GetNeutrino().CCNC(),mc.GetNeutrino().Mode()) << std::endl;
    std::cout << "-----------> Particles in the Stack = " << mc.NParticles() << std::endl;
    std::cout << std::setiosflags(std::ios::left) 
	      << std::setw(20) << "PARTICLE"
	      << std::setiosflags(std::ios::left) 
	      << std::setw(32) << "STATUS"
	      << std::setw(18) << "E (GeV)"
	      << std::setw(18) << "m (GeV/c2)"
	      << std::setw(18) << "Ek (GeV)"
	      << std::endl << std::endl;
    
    const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
    
    // Loop over the particle stack for this event 
    for(int i = 0; i < mc.NParticles(); ++i){
      simb::MCParticle part(mc.GetParticle(i));
      std::string name;
      if (part.PdgCode() == 18040)
	name = "Ar40 18040";
      else if (part.PdgCode() != -99999 )
	{
	  name = databasePDG->GetParticle(part.PdgCode())->GetName(); 
	}
      
      int code = part.StatusCode();
      std::string status = ParticleStatus(code);
      double mass = part.Mass();
      double energy = part.E(); 
      double Ek = (energy-mass); // Kinetic Energy (GeV)
      
      std::cout << std::setiosflags(std::ios::left) << std::setw(20) << name
		<< std::setiosflags(std::ios::left) << std::setw(32) <<status
		<< std::setw(18)<< energy
		<< std::setw(18)<< mass
		<< std::setw(18)<< Ek <<std::endl;  
    }


    if(mc.GetNeutrino().CCNC() == simb::kCC){
  
      ///look for the outgoing lepton in the particle stack
      ///just interested in the first one
      for(int i = 0; i < mc.NParticles(); ++i){
	simb::MCParticle part(mc.GetParticle(i));
	if(abs(part.PdgCode()) == 11){
	  fEMomentum->Fill(part.P());
	  fEDCosX->Fill(part.Px()/part.P());
	  fEDCosY->Fill(part.Py()/part.P());
	  fEDCosZ->Fill(part.Pz()/part.P());
	  fECons->Fill(nu.E() - part.E());
	  break;
	}
	else if(abs(part.PdgCode()) == 13){
	  fMuMomentum->Fill(part.P());
	  fMuDCosX->Fill(part.Px()/part.P());
	  fMuDCosY->Fill(part.Py()/part.P());
	  fMuDCosZ->Fill(part.Pz()/part.P());
	  fECons->Fill(nu.E() - part.E());
	  break;
	}
      }// end loop over particles
    }//end if CC interaction

    return;
  }

}

namespace evgen{

  DEFINE_ART_MODULE(NUANCEGen)

}

#endif 
