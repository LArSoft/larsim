////////////////////////////////////////////////////////////////////////
//
//
// NDK neutrino event generator
//
// echurch@fnal.gov
//
////////////////////////////////////////////////////////////////////////
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

#include "CLHEP/Random/RandFlat.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"


// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nutools/EventGeneratorBase/evgenbase.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"

#include "art/Framework/Core/EDProducer.h"

#include "TStopwatch.h"

class TH1F;
class TH2F;

namespace simb { class MCTruth;     }


namespace evgen {
  /// A module to check the results from the Monte Carlo generator
  class NDKGen : public art::EDProducer {
  public:
    explicit NDKGen(fhicl::ParameterSet const &pset);
    virtual ~NDKGen();                        

  private:

    void produce(art::Event& evt) override;
    void beginJob() override;
    void beginRun(art::Run& run) override;
    void endJob() override;

        std::string ParticleStatus(int StatusCode);
        std::string ReactionChannel(int ccnc,int mode);
    
    void FillHistograms(simb::MCTruth const& mc);

	std::string         fNdkFile;
    std::ifstream       fEventFile;
	TStopwatch          fStopwatch;      ///keep track of how long it takes to run the job
	
	std::string fNDKModuleLabel;
    CLHEP::HepRandomEngine& fEngine; ///< art-managed random-number engine
	
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
	
	TH1F* fECons;      ///< histogram to determine if energy is conserved in the event
	
  };
} // namespace

namespace evgen{

  //____________________________________________________________________________
  NDKGen::NDKGen(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fNdkFile{pset.get<std::string>("NdkFile")}
    , fEventFile{fNdkFile}
      // create a default random engine; obtain the random seed from NuRandomService,
      // unless overridden in configuration with key "Seed"
    , fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, pset, "Seed"))
  {
    fStopwatch.Start();

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();

    if(!fEventFile.good()) {
      throw cet::exception("NDKGen")
        << "Could not open file: " << fNdkFile << '\n';
    }

  }

  //____________________________________________________________________________
  NDKGen::~NDKGen()
  {  
    fStopwatch.Stop();
  }

//___________________________________________________________________________
  
  void NDKGen::beginJob()
  {
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
  void NDKGen::beginRun(art::Run& run)
  {
    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
  }

  //____________________________________________________________________________
  void NDKGen::endJob()
  {
    fEventFile.close();
  }

  //____________________________________________________________________________
  void NDKGen::produce(art::Event& evt)
  {

    std::cout << std::endl;
    std::cout<<"------------------------------------------------------------------------------"<<std::endl;
    //std::cout << "run    : " << evt.Header().Run() << std::endl;
    //std::cout << "subrun : " << evt.Header().Subrun() << std::endl;
    //std::cout << "event  : " << evt.Header().Event() << std::endl;
    std::cout << "event  : " << evt.id().event() << std::endl;  

    // TODO: fill more quantities out, as below.
    /*
    double X; // vertex position from Ndk
    double Y; // vertex position from Ndk
    double Z; // vertex position from Ndk
    double PDGCODE = -9999.; 
    double CHANNEL = -9999.;
    int channel = -9999;
    double energy = 0.; // in MeV from Ndk
    double cosx = 0.;
    double cosy = 0.;
    double cosz = 0.;

    int partnumber = 0;
    */

    std::string name, k, dollar;
     

    // event dump format on file output by the two commands ....
    // gevgen_ndcy  -g 1000180400  -m 8 -n 400 -o ndk
    // gevdump -f ndk.1000.ghep.root > ndk.out
    std::string Name;
    int Idx, Ist, PDG, Mother1, Mother2, Daughter1 ,Daughter2;
    double Px, Py, Pz, E, m ;
    std::string p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14;


    int trackid = -1; // set track id to -i as these are all primary particles and have id <= 0
    std::string primary("primary");
    int FirstMother = -1;
    double Mass = -9999;
    int Status = -9999;
    
    double P; // momentum of MCParticle IN GeV/c

    // TODO: Could perhaps imagine using these in NDk. 
    /*
    int targetnucleusPdg = -9999;
    int hitquarkPdg = -9999;
    double Q2 = -9999;
    */
    TLorentzVector Neutrino;
    TLorentzVector Lepton;
    TLorentzVector Target;
    TLorentzVector q;
    TLorentzVector Hadron4mom;

    // TODO: Could perhaps imagine using these in NDk. 
    /*
    int Tpdg = 0;  // for target 
    double Tmass = 0;
    int Tstatus = 11;
    double Tcosx, Tcosy, Tcosz, Tenergy;
    */

    TLorentzVector Tpos;
    
        
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
    simb::MCTruth truth;

    art::ServiceHandle<geo::Geometry> geo;
    CLHEP::RandFlat flat(fEngine);

    double const fvCut{5.0}; // force vtx to be this far from any wall.

    // Find boundary of active volume
    double minx = 1e9;
    double maxx = -1e9;
    double miny = 1e9;
    double maxy = -1e9;
    double minz = 1e9;
    double maxz = -1e9;
    for (size_t i = 0; i<geo->NTPC(); ++i)
    {
      double local[3] = {0.,0.,0.};
      double world[3] = {0.,0.,0.};
      const geo::TPCGeo &tpc = geo->TPC(i);
      tpc.LocalToWorld(local,world);
      if (minx>world[0]-geo->DetHalfWidth(i)) 
        minx = world[0]-geo->DetHalfWidth(i);
      if (maxx<world[0]+geo->DetHalfWidth(i))
        maxx = world[0]+geo->DetHalfWidth(i);
      if (miny>world[1]-geo->DetHalfHeight(i))
        miny = world[1]-geo->DetHalfHeight(i);
      if (maxy<world[1]+geo->DetHalfHeight(i))
        maxy = world[1]+geo->DetHalfHeight(i);
      if (minz>world[2]-geo->DetLength(i)/2.)
        minz = world[2]-geo->DetLength(i)/2.;
      if (maxz<world[2]+geo->DetLength(i)/2.)
        maxz = world[2]+geo->DetLength(i)/2.;
    }

    // Assign vertice position
    double X0 = 0.0 + flat.fire( minx+fvCut , maxx-fvCut );
    double Y0 = 0.0 + flat.fire( miny+fvCut , maxy-fvCut );
    double Z0 = 0.0 + flat.fire( minz+fvCut , maxz-fvCut );

    std::cout << "NDKGen_module: X, Y, Z of vtx: " << X0 << ", "<< Y0 << ", "<< Z0 << std::endl;
    
    int GenieEvt = -999;
    
    if(!fEventFile.good())
      std::cout << "NdkFile: Problem reading Ndk file" << std::endl; 
    
    while(getline(fEventFile,k)){

      if (k.find("** Event:")!= std::string::npos) {
        std::istringstream in;
        in.clear();
        in.str(k);
        std::string dummy;   
        in>> dummy>> dummy>> dummy >> dummy>> dummy>> dummy>> dummy >> dummy>> dummy>> dummy >> GenieEvt;
        std::cout<<"Genie Evt "<< GenieEvt <<" art evt "<<evt.id().event()<<"\n";
      }

      if (GenieEvt+1 != static_cast<int>(evt.id().event()))
        continue;
      else {

        if (!k.compare(0,25,"GENIE Interaction Summary")) // testing for new event.
          break;
        if (k.compare(0,1,"|") || k.compare(1,2,"  ")) continue; // uninteresting line if it doesn't start with "|" and if second and third characters aren't spaces.
        if (k.find("Fin-Init") != std::string::npos) continue; // Meh.
        if (k.find("Ar") != std::string::npos) continue; // Meh.
        if (k.find("Cl") != std::string::npos) continue; // ignore chlorine nucleus in nnbar events
        if (k.find("HadrBlob") != std::string::npos) continue; // Meh.
        if (k.find("NucBindE") != std::string::npos) continue; // Meh. atmo
        if (k.find("FLAGS") != std::string::npos) break; // Event end. gevgen_ndcy
        if (k.find("Vertex") != std::string::npos) break; // Event end. atmo

        //      if (!k.compare(26,1,"3") || !k.compare(26,1,"1")) ; // New event or stable particles.
        if (!k.compare(26,1,"1"))  // New event or stable particles.
        {

      	  std::istringstream in;
      	  in.clear();
      	  in.str(k);

      	  in>>p1>> Idx >>p2>> Name >>p3>> Ist >>p4>> PDG >>p5>>Mother1 >> p6 >> Mother2 >>p7>> Daughter1 >>p8>> Daughter2 >>p9>>Px>>p10>>Py>>p11>>Pz>>p12>>E>>p13>> m>>p14;
          //std::cout<<std::setprecision(9)<<dollar<<"  "<<name<<"  "<<PDGCODE<<"  "<<energy<<"  "<<cosx<<" "<<cosy<<"  "<<cosz<<"  "<<partnumber<<std::endl;
      	  if (Ist!=1) continue;

      	  std::cout << "PDG = " << PDG << std::endl;	
      	
      	  const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
      	  const TParticlePDG* definition = databasePDG->GetParticle(PDG);
      	  Mass = definition->Mass(); // GeV
      	  if (E-Mass < 0.001) continue; // KE is too low.
      	
      	  //	  if(partnumber == -1)
      	  Status = 1;
      	  
      	  simb::MCParticle mcpart(trackid,
      		  		PDG,
      			  	primary,		    
      		  		FirstMother,
      			  	Mass,
      		  		Status
      			  	);
      		
      	  P = std::sqrt(pow(E,2.) - pow(Mass,2.)); // GeV/c
      	  std::cout << "Momentum = " << P << std::endl;
      	
      	  TLorentzVector pos(X0, Y0, Z0, 0);

      	  Tpos = pos; // for target
      	
      	  TLorentzVector mom(Px,Py,Pz, E);
      	
      	  mcpart.AddTrajectoryPoint(pos,mom);
      	  truth.Add(mcpart);
      	
      	
        }// loop over particles in an event
        truth.SetOrigin(simb::kUnknown);
            
        //if (!k.compare(1,1,"FLAGS")) // end of event
        //  break;  

      }
    } // end while loop
    
    /////////////////////////////////
    std::cout << "NDKGen.cxx: Putting " << truth.NParticles() << " tracks on stack." << std::endl; 
    truthcol->push_back(truth);
    //FillHistograms(truth);  
    evt.put(std::move(truthcol));
      
    return;
  }
  
//   //......................................................................
  std::string NDKGen::ParticleStatus(int StatusCode)
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
  std::string NDKGen::ReactionChannel(int ccnc,int mode)
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
  void NDKGen::FillHistograms(simb::MCTruth const& mc)
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


//     MF_LOG_DEBUG("GENIEInteractionInformation") 
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
// 	MF_LOG_DEBUG("GENIEFinalState")
// 	  << std::setiosflags(std::ios::left) << std::setw(20) << name
// 	  << std::setiosflags(std::ios::left) << std::setw(32) <<status
// 	  << std::setw(18)<< energy
// 	  << std::setw(18)<< mass
// 	  << std::setw(18)<< Ek <<std::endl;
//       else 
// 	MF_LOG_DEBUG("GENIEFinalState") 
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

  DEFINE_ART_MODULE(NDKGen)

} // namespace
