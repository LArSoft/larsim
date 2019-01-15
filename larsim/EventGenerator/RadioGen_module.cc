////////////////////////////////////////////////////////////////////////
/// \file  RadioGen_module.cc
/// \brief Generator for radiological decays
///
/// Module designed to produce a set list of particles for MC to model radiological decays
///
/// \author  trj@fnal.gov
//           Rn222 generation feature added by gleb.sinev@duke.edu
//           (based on a generator by jason.stock@mines.sdsmt.edu)
//           JStock. Added preliminary changes to get ready for Ar42 implimentation. This includes allowing for multiple particles from the same decay.
//           Ar42 generation added by JStock (jason.stock@mines.sdsmt.edu).
//             Ar42 is designed to handle 5 separate different
//             beta decay modes, each with it's own chains of
//             possible dexcitation gammas. Because of the
//             high energies, and the relatively high rate
//             expected in the DUNE FD, these chains are
//             completely simulated instead of relying on the
//             existing machinery in this module. To make the
//             treatment of multiple decay products and
//             dexcitation chains generally available would
//             require a significant redesign of the module
//             and possibly of the .root files data structure
//             for each radiological.
//
//           Dec 01, 2017 JStock
//             Adding the ability to make 8.997 MeV Gammas for
//             Ni59 Calibration sources. This is another "hacky"
//             fix to something that really deserves a more
//             elegant and comprehensive solution.
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_RADIOLOGICAL
#define EVGEN_RADIOLOGICAL

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <regex>
#include <cmath>
#include <memory>
#include <iterator>
#include <vector>
#include <sys/stat.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoNode.h>

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
#include "cetlib/search_path.h"
#include "cetlib/exempt_ptr.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"

// root includes

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGenPhaseSpace.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraph.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

namespace simb { class MCTruth; }

namespace evgen {

  /// Module to generate particles created by radiological decay, patterend off of SingleGen
  /// Currently it generates only in rectangular prisms oriented along the x,y,z axes



  class RadioGen : public art::EDProducer {

  public:
    explicit RadioGen(fhicl::ParameterSet const& pset);
    virtual ~RadioGen();

    // This is called for each event.
    void produce(art::Event& evt);
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    typedef int    ti_PDGID;  // These typedefs may look odd, and unecessary. I chose to use them to make the tuples I use later more readable. ti, type integer :JStock
    typedef double td_Mass;   // These typedefs may look odd, and unecessary. I chose to use them to make the tuples I use later more readable. td, type double  :JStock

    void SampleOne(unsigned int   i,
       simb::MCTruth &mct);

    TLorentzVector dirCalc(double p, double m);

    void readfile(std::string nuclide, std::string filename);
    void samplespectrum(std::string nuclide, int &itype, double &t, double &m, double &p);

    void Ar42Gamma2(std::vector<std::tuple<ti_PDGID, td_Mass, TLorentzVector>>& v_prods);
    void Ar42Gamma3(std::vector<std::tuple<ti_PDGID, td_Mass, TLorentzVector>>& v_prods);
    void Ar42Gamma4(std::vector<std::tuple<ti_PDGID, td_Mass, TLorentzVector>>& v_prods);
    void Ar42Gamma5(std::vector<std::tuple<ti_PDGID, td_Mass, TLorentzVector>>& v_prods);

    // recoded so as to use the LArSoft-managed random number generator
    double samplefromth1d(TH1D *hist);

    // itype = pdg code:  1000020040: alpha.  itype=11: beta. -11: positron,  itype=22: gamma.  -1: error
    // t is the kinetic energy in GeV  (= E-m)
    // m = mass in GeV
    // p = momentum in GeV

    //double betaphasespace(double mass, double q); // older parameterization.

    // the generator randomly samples points in a rectangular prism of space and time, and only selects those points in
    // volumes with materials that match the regexes in fMaterial.  One can use wildcards * and ? for broader matches.

    std::vector<std::string> fNuclide;   ///< List of nuclides to simulate.  Example:  "39Ar".
    std::vector<std::string> fMaterial;  ///< List of regexes of materials in which to generate the decays.  Example: "LAr"
    std::vector<double> fBq;             ///< Radioactivity in Becquerels (decay per sec) per cubic cm.
    std::vector<double> fT0;             ///< Beginning of time window to simulate in ns
    std::vector<double> fT1;             ///< End of time window to simulate in ns
    std::vector<double> fX0;             ///< Bottom corner x position (cm) in world coordinates
    std::vector<double> fY0;             ///< Bottom corner y position (cm) in world coordinates
    std::vector<double> fZ0;             ///< Bottom corner z position (cm) in world coordinates
    std::vector<double> fX1;             ///< Top corner x position (cm) in world coordinates
    std::vector<double> fY1;             ///< Top corner y position (cm) in world coordinates
    std::vector<double> fZ1;             ///< Top corner z position (cm) in world coordinates
    bool                fIsFirstSignalSpecial;
    int trackidcounter;                  ///< Serial number for the MC track ID


    // leftovers from the phase space generator
    // const double gevperamu = 0.931494061;
    // TGenPhaseSpace rg;  // put this here so we don't constantly construct and destruct it

    const double m_e = 0.000510998928;  // mass of electron in GeV
    const double m_alpha = 3.727379240; // mass of an alpha particle in GeV
    const double m_neutron = 0.9395654133; // mass of a neutron in GeV

    std::vector<std::string> spectrumname;
    std::vector<TH1D*> alphaspectrum;
    std::vector<double> alphaintegral;
    std::vector<TH1D*> betaspectrum;
    std::vector<double> betaintegral;
    std::vector<TH1D*> gammaspectrum;
    std::vector<double> gammaintegral;
    std::vector<TH1D*> neutronspectrum;
    std::vector<double> neutronintegral;
    cet::exempt_ptr<CLHEP::HepRandomEngine> fEngine{nullptr};
  };
}

namespace evgen{

  //____________________________________________________________________________
  RadioGen::RadioGen(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
  {

    this->reconfigure(pset);

    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<rndm::NuRandomService>()
      ->createEngine(*this, pset, "Seed");
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    auto& engine = rng->getEngine(art::ScheduleID::first(),
                                  pset.get<std::string>("module_label"));
    fEngine = cet::make_exempt_ptr(&engine);

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();

  }

  //____________________________________________________________________________
  RadioGen::~RadioGen()
  {
  }

  //____________________________________________________________________________
  void RadioGen::reconfigure(fhicl::ParameterSet const& p)
  {
    // do not put seed in reconfigure because we don't want to reset
    // the seed midstream -- same as SingleGen

    fNuclide       = p.get< std::vector<std::string>>("Nuclide");
    fMaterial      = p.get< std::vector<std::string>>("Material");
    fBq            = p.get< std::vector<double> >("BqPercc");
    fT0            = p.get< std::vector<double> >("T0");
    fT1            = p.get< std::vector<double> >("T1");
    fX0            = p.get< std::vector<double> >("X0");
    fY0            = p.get< std::vector<double> >("Y0");
    fZ0            = p.get< std::vector<double> >("Z0");
    fX1            = p.get< std::vector<double> >("X1");
    fY1            = p.get< std::vector<double> >("Y1");
    fZ1            = p.get< std::vector<double> >("Z1");
    fIsFirstSignalSpecial      = p.get< bool >("IsFirstSignalSpecial", false);

    // check for consistency of vector sizes

    unsigned int nsize = fNuclide.size();
    if (  fMaterial.size() != nsize ) throw cet::exception("RadioGen") << "Different size Material vector and Nuclide vector\n";
    if (  fBq.size() != nsize ) throw cet::exception("RadioGen") << "Different size Bq vector and Nuclide vector\n";
    if (  fT0.size() != nsize ) throw cet::exception("RadioGen") << "Different size T0 vector and Nuclide vector\n";
    if (  fT1.size() != nsize ) throw cet::exception("RadioGen") << "Different size T1 vector and Nuclide vector\n";
    if (  fX0.size() != nsize ) throw cet::exception("RadioGen") << "Different size X0 vector and Nuclide vector\n";
    if (  fY0.size() != nsize ) throw cet::exception("RadioGen") << "Different size Y0 vector and Nuclide vector\n";
    if (  fZ0.size() != nsize ) throw cet::exception("RadioGen") << "Different size Z0 vector and Nuclide vector\n";
    if (  fX1.size() != nsize ) throw cet::exception("RadioGen") << "Different size X1 vector and Nuclide vector\n";
    if (  fY1.size() != nsize ) throw cet::exception("RadioGen") << "Different size Y1 vector and Nuclide vector\n";
    if (  fZ1.size() != nsize ) throw cet::exception("RadioGen") << "Different size Z1 vector and Nuclide vector\n";

    for(std::string & nuclideName : fNuclide){
      if(nuclideName=="39Ar"      ){readfile("39Ar","Argon_39.root")    ;}
      else if(nuclideName=="60Co" ){readfile("60Co","Cobalt_60.root")   ;}
      else if(nuclideName=="85Kr" ){readfile("85Kr","Krypton_85.root")  ;}
      else if(nuclideName=="40K"  ){readfile("40K","Potassium_40.root") ;}
      else if(nuclideName=="232Th"){readfile("232Th","Thorium_232.root");}
      else if(nuclideName=="238U" ){readfile("238U","Uranium_238.root") ;}
      else if(nuclideName=="222Rn"){continue;} //Rn222 is handeled separately later
      else if(nuclideName=="59Ni"){continue;} //Rn222 is handeled separately later
      else if(nuclideName=="42Ar" ){
        readfile("42Ar_1", "Argon_42_1.root"); //Each possible beta decay mode of Ar42 is given it's own .root file for now.
        readfile("42Ar_2", "Argon_42_2.root"); //This allows us to know which decay chain to follow for the dexcitation gammas.
        readfile("42Ar_3", "Argon_42_3.root"); //The dexcitation gammas are not included in the root files as we want to
        readfile("42Ar_4", "Argon_42_4.root"); //probabilistically simulate the correct coincident gammas, which we cannot guarantee
        readfile("42Ar_5", "Argon_42_5.root"); //by sampling a histogram.
        continue;
      } //Ar42  is handeled separately later
      else{
        std::string searchName = nuclideName;
        searchName+=".root";
        readfile(nuclideName,searchName);
      }
    }

    return;
  }

  //____________________________________________________________________________
  void RadioGen::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using

    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
    run.put(std::move(runcol));

    return;
  }

  //____________________________________________________________________________
  void RadioGen::produce(art::Event& evt)
  {

    ///unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);

    trackidcounter = -1;
    for (unsigned int i=0; i<fNuclide.size(); ++i) {
      SampleOne(i,truth);
    }//end loop over nuclides

    MF_LOG_DEBUG("RadioGen") << truth;
    truthcol->push_back(truth);
    evt.put(std::move(truthcol));
    return;
  }

  //____________________________________________________________________________
  // Generate radioactive decays per nuclide per volume according to the FCL parameters

  void RadioGen::SampleOne(unsigned int i, simb::MCTruth &mct){

    art::ServiceHandle<geo::Geometry> geo;
    TGeoManager *geomanager = geo->ROOTGeoManager();

    // get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::RandFlat     flat(*fEngine);
    CLHEP::RandPoisson  poisson(*fEngine);

    // figure out how many decays to generate, assuming that the entire prism consists of the radioactive material.
    // we will skip over decays in other materials later.

    double rate = fabs( fBq[i] * (fT1[i] - fT0[i]) * (fX1[i] - fX0[i]) * (fY1[i] - fY0[i]) * (fZ1[i] - fZ0[i]) ) / 1.0E9;
    long ndecays = poisson.shoot(rate);

    for (unsigned int idecay=0; idecay<ndecays; idecay++)
    {
      // generate just one particle at a time.  Need a little recoding if a radioactive
      // decay generates multiple daughters that need simulation
      // uniformly distributed in position and time
      //
      // JStock: Leaving this as a single position for the decay products. For now I will assume they all come from the same spot.
      TLorentzVector pos( fX0[i] + flat.fire()*(fX1[i] - fX0[i]),
          fY0[i] + flat.fire()*(fY1[i] - fY0[i]),
          fZ0[i] + flat.fire()*(fZ1[i] - fZ0[i]),
          (idecay==0 && fIsFirstSignalSpecial) ? 0 : ( fT0[i] + flat.fire()*(fT1[i] - fT0[i]) ) );
          //fT0[i] + flat.fire()*(fT1[i] - fT0[i]) );

      // discard decays that are not in the proper material
      std::string volmaterial = geomanager->FindNode(pos.X(),pos.Y(),pos.Z())->GetMedium()->GetMaterial()->GetName();
      //mf::LogDebug("RadioGen") << "Position: " << pos.X() << " " << pos.Y() << " " << pos.Z() << " " << pos.T() << " " << volmaterial << " " << fMaterial[i] << std::endl;
      if ( ! std::regex_match(volmaterial, std::regex(fMaterial[i])) ) continue;
      //mf::LogDebug("RadioGen") << "Decay accepted" << std::endl;

      //Moved pdgid into the next statement, so that it is localized.
      // electron=11, photon=22, alpha = 1000020040, neutron = 2112

      //JStock: Allow us to have different particles from the same decay. This requires multiple momenta.
      std::vector<std::tuple<ti_PDGID, td_Mass, TLorentzVector>> v_prods; //(First is for PDGID, second is mass, third is Momentum)

      if (fNuclide[i] == "222Rn")          // Treat 222Rn separately
      {
        double p=0; double t=0.00548952; td_Mass m=m_alpha; ti_PDGID pdgid=1000020040; //td_Mass = double. ti_PDGID = int;
        double energy = t + m;
        double p2     = energy*energy - m*m;
        if (p2 > 0) p = TMath::Sqrt(p2);
        else        p = 0;
        //Make TLorentzVector and push it to the back of v_prods.
        std::tuple<ti_PDGID, td_Mass, TLorentzVector> partMassMom = std::make_tuple(pdgid, m, dirCalc(p,m));
        v_prods.push_back(partMassMom);
      }//End special case RN222
      else if(fNuclide[i] == "59Ni"){ //Treat 59Ni Calibration Source separately (as I haven't made a spectrum for it, and ultimately it should be handeled with multiple particle outputs.
        double p=0.008997; td_Mass m=0; ti_PDGID pdgid=22; // td_Mas=double. ti_PDFID=int. Assigning p directly, as t=p for gammas.
        std::tuple<ti_PDGID, td_Mass, TLorentzVector> partMassMom = std::make_tuple(pdgid, m, dirCalc(p,m));
        v_prods.push_back(partMassMom);
      }//end special case Ni59 calibration source
      else if(fNuclide[i] == "42Ar"){   // Spot for special treatment of Ar42.
        double p=0; double t=0; td_Mass m = 0; ti_PDGID pdgid=0; //td_Mass = double. ti_PDGID = int;
        double bSelect = flat.fire();   //Make this a random number from 0 to 1.
        if(bSelect<0.819){              //beta channel 1. No Gamma. beta Q value 3525.22 keV
          samplespectrum("42Ar_1", pdgid, t, m, p);
          std::tuple<ti_PDGID, td_Mass, TLorentzVector> partMassMom = std::make_tuple(pdgid, m, dirCalc(p,m));
          v_prods.push_back(partMassMom);
          //No gamma here.
        }else if(bSelect<0.9954){       //beta channel 2. 1 Gamma (1524.6 keV). beta Q value 2000.62
          samplespectrum("42Ar_2", pdgid, t, m, p);
          std::tuple<ti_PDGID, td_Mass, TLorentzVector> partMassMom = std::make_tuple(pdgid, m, dirCalc(p,m));
          v_prods.push_back(partMassMom);
          Ar42Gamma2(v_prods);
        }else if(bSelect<0.9988){       //beta channel 3. 1 Gamma Channel. 312.6 keV + gamma 2. beta Q value 1688.02 keV
          samplespectrum("42Ar_3", pdgid, t, m, p);
          std::tuple<ti_PDGID, td_Mass, TLorentzVector> partMassMom = std::make_tuple(pdgid, m, dirCalc(p,m));
          v_prods.push_back(partMassMom);
          Ar42Gamma3(v_prods);
        }else if(bSelect<0.9993){       //beta channel 4. 2 Gamma Channels. Either 899.7 keV (i 0.052) + gamma 2 or 2424.3 keV (i 0.020). beta Q value 1100.92 keV
          samplespectrum("42Ar_4", pdgid, t, m, p);
          std::tuple<ti_PDGID, td_Mass, TLorentzVector> partMassMom = std::make_tuple(pdgid, m, dirCalc(p,m));
          v_prods.push_back(partMassMom);
          Ar42Gamma4(v_prods);
        }else{                          //beta channel 5. 3 gamma channels. 692.0 keV + 1228.0 keV + Gamma 2 (i 0.0033) ||OR|| 1021.2 keV + gamma 4 (i 0.0201) ||OR|| 1920.8 keV + gamma 2 (i 0.041). beta Q value 79.82 keV
          samplespectrum("42Ar_5", pdgid, t, m, p);
          std::tuple<ti_PDGID, td_Mass, TLorentzVector> partMassMom = std::make_tuple(pdgid, m, dirCalc(p,m));
          v_prods.push_back(partMassMom);
          Ar42Gamma5(v_prods);
        }
        //Add beta.
        //Call gamma function for beta mode.
      }
      else{ //General Case.
        double p=0; double t=0; td_Mass m = 0; ti_PDGID pdgid=0; //td_Mass = double. ti_PDGID = int;
        samplespectrum(fNuclide[i],pdgid,t,m,p);
        std::tuple<ti_PDGID, td_Mass, TLorentzVector> partMassMom = std::make_tuple(pdgid, m, dirCalc(p,m));
        v_prods.push_back(partMassMom);
      }//end else (not RN or other special case

      //JStock: Modify this to now loop over the v_prods.
      for(auto prodEntry : v_prods){
        // set track id to a negative serial number as these are all primary particles and have id <= 0
        int trackid = trackidcounter;
        ti_PDGID pdgid = std::get<0>(prodEntry);
        td_Mass  m = std::get<1>(prodEntry);
        TLorentzVector pvec = std::get<2>(prodEntry);
        trackidcounter--;
        std::string primary("primary");

        // alpha particles need a little help since they're not in the TDatabasePDG table
        // // so don't rely so heavily on default arguments to the MCParticle constructor
        if (pdgid == 1000020040){
          simb::MCParticle part(trackid, pdgid, primary,-1,m,1);
          part.AddTrajectoryPoint(pos, pvec);
          mct.Add(part);
        }// end "If alpha"
        else{
          simb::MCParticle part(trackid, pdgid, primary);
          part.AddTrajectoryPoint(pos, pvec);
          mct.Add(part);
        }// end All standard cases.
      }//End Loop over all particles produces in this single decay.
    }
  }

  //Calculate an arbitrary direction with a given magnitude p
  TLorentzVector RadioGen::dirCalc(double p, double m){
      CLHEP::RandFlat  flat(*fEngine);
      // isotropic production angle for the decay product
      double costheta = (2.0*flat.fire() - 1.0);
      if (costheta < -1.0) costheta = -1.0;
      if (costheta > 1.0) costheta = 1.0;
      double sintheta = sqrt(1.0-costheta*costheta);
      double phi = 2.0*M_PI*flat.fire();
        TLorentzVector pvec(p*sintheta*std::cos(phi),
          p*sintheta*std::sin(phi),
          p*costheta,
          std::sqrt(p*p+m*m));
      return pvec;
  }


  // only reads those files that are on the fNuclide list.  Copy information from the TGraphs to TH1D's

  void RadioGen::readfile(std::string nuclide, std::string filename)
  {
    int ifound = 0;
    for (size_t i=0; i<fNuclide.size(); i++)
    {
      if (fNuclide[i] == nuclide){ //This check makes sure that the nuclide we are searching for is in fact in our fNuclide list. Ar42 handeled separately.
        ifound = 1;
        break;
      } //End If nuclide is in our list. Next is the special case of Ar42
      else if ( (std::regex_match(nuclide, std::regex( "42Ar.*" )) ) && fNuclide[i]=="42Ar" ){
        ifound = 1;
        break;
      }
    }
    if (ifound == 0) return;

    Bool_t addStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE); // cloned histograms go in memory, and aren't deleted when files are closed.
    // be sure to restore this state when we're out of the routine.


    spectrumname.push_back(nuclide);
    cet::search_path sp("FW_SEARCH_PATH");
    std::string fn2 = "Radionuclides/";
    fn2 += filename;
    std::string fullname;
    sp.find_file(fn2, fullname);
    struct stat sb;
    if (fullname.empty() || stat(fullname.c_str(), &sb)!=0)
      throw cet::exception("RadioGen") << "Input spectrum file "
        << fn2
        << " not found in FW_SEARCH_PATH!\n";

    TFile f(fullname.c_str(),"READ");
    TGraph *alphagraph = (TGraph*) f.Get("Alphas");
    TGraph *betagraph = (TGraph*) f.Get("Betas");
    TGraph *gammagraph = (TGraph*) f.Get("Gammas");
    TGraph *neutrongraph = (TGraph*) f.Get("Neutrons");

    if (alphagraph)
    {
      int np = alphagraph->GetN();
      double *y = alphagraph->GetY();
      std::string name;
      name = "RadioGen_";
      name += nuclide;
      name += "_Alpha";
      TH1D *alphahist = (TH1D*) new TH1D(name.c_str(),"Alpha Spectrum",np,0,np);
      for (int i=0; i<np; i++)
      {
        alphahist->SetBinContent(i+1,y[i]);
        alphahist->SetBinError(i+1,0);
      }
      alphaspectrum.push_back(alphahist);
      alphaintegral.push_back(alphahist->Integral());
    }
    else
    {
      alphaspectrum.push_back(0);
      alphaintegral.push_back(0);
    }


    if (betagraph)
    {
      int np = betagraph->GetN();

      double *y = betagraph->GetY();
      std::string name;
      name = "RadioGen_";
      name += nuclide;
      name += "_Beta";
      TH1D *betahist = (TH1D*) new TH1D(name.c_str(),"Beta Spectrum",np,0,np);

      for (int i=0; i<np; i++)
      {
        betahist->SetBinContent(i+1,y[i]);
        betahist->SetBinError(i+1,0);
      }
      betaspectrum.push_back(betahist);
      betaintegral.push_back(betahist->Integral());
    }
    else
    {
      betaspectrum.push_back(0);
      betaintegral.push_back(0);
    }

    if (gammagraph)
    {
      int np = gammagraph->GetN();
      double *y = gammagraph->GetY();
      std::string name;
      name = "RadioGen_";
      name += nuclide;
      name += "_Gamma";
      TH1D *gammahist = (TH1D*) new TH1D(name.c_str(),"Gamma Spectrum",np,0,np);
      for (int i=0; i<np; i++)
      {
        gammahist->SetBinContent(i+1,y[i]);
        gammahist->SetBinError(i+1,0);
      }
      gammaspectrum.push_back(gammahist);
      gammaintegral.push_back(gammahist->Integral());
    }
    else
    {
      gammaspectrum.push_back(0);
      gammaintegral.push_back(0);
    }

    if (neutrongraph)
    {
      int np = neutrongraph->GetN();
      double *y = neutrongraph->GetY();
      std::string name;
      name = "RadioGen_";
      name += nuclide;
      name += "_Neutron";
      TH1D *neutronhist = (TH1D*) new TH1D(name.c_str(),"Neutron Spectrum",np,0,np);
      for (int i=0; i<np; i++)
      {
        neutronhist->SetBinContent(i+1,y[i]);
        neutronhist->SetBinError(i+1,0);
      }
      neutronspectrum.push_back(neutronhist);
      neutronintegral.push_back(neutronhist->Integral());
    }
    else
    {
      neutronspectrum.push_back(0);
      neutronintegral.push_back(0);
    }

    f.Close();
    TH1::AddDirectory(addStatus);

    double total = alphaintegral.back() + betaintegral.back() + gammaintegral.back() + neutronintegral.back();
    if (total>0)
    {
      alphaintegral.back() /= total;
      betaintegral.back() /= total;
      gammaintegral.back() /= total;
      neutronintegral.back() /= total;
    }
  }


  void RadioGen::samplespectrum(std::string nuclide, int &itype, double &t, double &m, double &p)
  {
    CLHEP::RandFlat  flat(*fEngine);

    int inuc = -1;
    for (size_t i=0; i<spectrumname.size(); i++)
    {
      if (nuclide == spectrumname[i])
      {
        inuc = i;
        break;
      }
    }
    if (inuc == -1)
    {
      t=0;  // throw an exception in the future
      itype = 0;
      throw cet::exception("RadioGen") << "Ununderstood nuclide:  " << nuclide << "\n";
    }

    double rtype = flat.fire();

    itype = -1;
    m = 0;
    p = 0;
    for (int itry=0;itry<10;itry++) // maybe a tiny normalization issue with a sum of 0.99999999999 or something, so try a few times.
    {
      if (rtype <= alphaintegral[inuc] && alphaspectrum[inuc] != 0)
      {
        itype = 1000020040; // alpha
        m = m_alpha;
        t = samplefromth1d(alphaspectrum[inuc])/1000000.0;
      }
      else if (rtype <= alphaintegral[inuc]+betaintegral[inuc] && betaspectrum[inuc] != 0)
      {
        itype = 11; // beta
        m = m_e;
        t = samplefromth1d(betaspectrum[inuc])/1000000.0;
      }
      else if ( rtype <= alphaintegral[inuc] + betaintegral[inuc] + gammaintegral[inuc] && gammaspectrum[inuc] != 0)
      {
        itype = 22; // gamma
        m = 0;
        t = samplefromth1d(gammaspectrum[inuc])/1000000.0;
      }
      else if( neutronspectrum[inuc] != 0)
      {
        itype = 2112;
        m     = m_neutron;
        t     = samplefromth1d(neutronspectrum[inuc])/1000000.0;
      }
      if (itype >= 0) break;
    }
    if (itype == -1)
    {
      throw cet::exception("RadioGen") << "Normalization problem with nuclide:  " << nuclide << "\n";
    }
    double e = t + m;
    p = e*e - m*m;
    if (p>=0)
    { p = TMath::Sqrt(p); }
    else
    { p=0; }
  }

  // this is just a copy of TH1::GetRandom that uses the art-managed CLHEP random number generator instead of gRandom
  // and a better handling of negative bin contents

  double RadioGen::samplefromth1d(TH1D *hist)
  {
    CLHEP::RandFlat  flat(*fEngine);

    int nbinsx = hist->GetNbinsX();
    std::vector<double> partialsum;
    partialsum.resize(nbinsx+1);
    partialsum[0] = 0;

    for (int i=1;i<=nbinsx;i++)
    {
      double hc = hist->GetBinContent(i);
      if ( hc < 0) throw cet::exception("RadioGen") << "Negative bin:  " << i << " " << hist->GetName() << "\n";
      partialsum[i] = partialsum[i-1] + hc;
    }
    double integral = partialsum[nbinsx];
    if (integral == 0) return 0;
    // normalize to unit sum
    for (int i=1;i<=nbinsx;i++) partialsum[i] /= integral;

    double r1 = flat.fire();
    int ibin = TMath::BinarySearch(nbinsx,&(partialsum[0]),r1);
    Double_t x = hist->GetBinLowEdge(ibin+1);
    if (r1 > partialsum[ibin]) x +=
      hist->GetBinWidth(ibin+1)*(r1-partialsum[ibin])/(partialsum[ibin+1] - partialsum[ibin]);
    return x;
  }


  //Ar42 uses BNL tables for K-42 from Aug 2017
  //beta channel 1. No Gamma. beta Q value 3525.22 keV
  //beta channel 2. 1 Gamma (1524.6 keV). beta Q value 2000.62
  //beta channel 3. 1 Gamma Channel. 312.6 keV + gamma 2. beta Q value 1688.02 keV
  //beta channel 4. 2 Gamma Channels. Either 899.7 keV (i 0.052) + gamma 2 or 2424.3 keV (i 0.020). beta Q value 1100.92 keV
  //beta channel 5. 3 gamma channels. 692.0 keV + 1228.0 keV + Gamma 2 (i 0.0033) ||OR|| 1021.2 keV + gamma 4 (i 0.0201) ||OR|| 1920.8 keV + gamma 2 (i 0.041). beta Q value 79.82 keV
  //No Ar42Gamma1 as beta channel 1 does not produce a dexcitation gamma.
  void RadioGen::Ar42Gamma2(std::vector<std::tuple<ti_PDGID, td_Mass, TLorentzVector>>& v_prods){
    ti_PDGID pdgid = 22; td_Mass m = 0.0; //we are writing gammas
    std::vector<double> vd_p = {.0015246};//Momentum in GeV
    for(auto p : vd_p){
      std::tuple<ti_PDGID, td_Mass, TLorentzVector> prods = std::make_tuple(pdgid, m, dirCalc(p,m));
      v_prods.push_back(prods);
    }
  }
  void RadioGen::Ar42Gamma3(std::vector<std::tuple<ti_PDGID, td_Mass, TLorentzVector>>& v_prods){
    ti_PDGID pdgid = 22; td_Mass m = 0.0; //we are writing gammas
    std::vector<double> vd_p = {.0003126};
    for(auto p : vd_p){
      std::tuple<ti_PDGID, td_Mass, TLorentzVector> prods = std::make_tuple(pdgid, m, dirCalc(p,m));
      v_prods.push_back(prods);
    }
    Ar42Gamma2(v_prods);
  }
  void RadioGen::Ar42Gamma4(std::vector<std::tuple<ti_PDGID, td_Mass, TLorentzVector>>& v_prods){
    CLHEP::RandFlat     flat(*fEngine);
    ti_PDGID pdgid = 22; td_Mass m = 0.0; //we are writing gammas
    double chan1 = (0.052 / (0.052+0.020) );
    if(flat.fire()<chan1){
      std::vector<double> vd_p = {.0008997};//Momentum in GeV
      for(auto p : vd_p){
        std::tuple<ti_PDGID, td_Mass, TLorentzVector> prods = std::make_tuple(pdgid, m, dirCalc(p,m));
        v_prods.push_back(prods);
      }
      Ar42Gamma2(v_prods);
    }else{
      std::vector<double> vd_p = {.0024243};//Momentum in GeV
      for(auto p : vd_p){
        std::tuple<ti_PDGID, td_Mass, TLorentzVector> prods = std::make_tuple(pdgid, m, dirCalc(p,m));
        v_prods.push_back(prods);
      }
    }
  }
  void RadioGen::Ar42Gamma5(std::vector<std::tuple<ti_PDGID, td_Mass, TLorentzVector>>& v_prods){
    CLHEP::RandFlat     flat(*fEngine);
    ti_PDGID pdgid = 22; td_Mass m = 0.0; //we are writing gammas
    double chan1 = ( 0.0033 / (0.0033 + 0.0201 + 0.041) ); double chan2 = ( 0.0201 / (0.0033 + 0.0201 + 0.041) );
    double chanPick = flat.fire();
    if(chanPick < chan1){
      std::vector<double> vd_p = {0.000692, 0.001228};//Momentum in GeV
      for(auto p : vd_p){
        std::tuple<ti_PDGID, td_Mass, TLorentzVector> prods = std::make_tuple(pdgid, m, dirCalc(p,m));
        v_prods.push_back(prods);
      }
      Ar42Gamma2(v_prods);
    }else if (chanPick<(chan1+chan2)){
      std::vector<double> vd_p = {0.0010212};//Momentum in GeV
      for(auto p : vd_p){
        std::tuple<ti_PDGID, td_Mass, TLorentzVector> prods = std::make_tuple(pdgid, m, dirCalc(p,m));
        v_prods.push_back(prods);
      }
      Ar42Gamma4(v_prods);
    }else{
      std::vector<double> vd_p = {0.0019208};//Momentum in GeV
      for(auto p : vd_p){
        std::tuple<ti_PDGID, td_Mass, TLorentzVector> prods = std::make_tuple(pdgid, m, dirCalc(p,m));
        v_prods.push_back(prods);
      }
      Ar42Gamma2(v_prods);
    }
  }

  // phase space generator for beta decay -- keep it as a comment in case we ever want to revive it

  // double RadioGen::betaphasespace(double mass, double q)
  //{
  //  art::ServiceHandle<art::RandomNumberGenerator> rng;
  //  CLHEP::HepRandomEngine &engine = rng->getEngine();
  //  CLHEP::RandFlat     flat(engine);
  //  double p = 0;
  //  double mi = mass+q+m_e;
  //  TLorentzVector p0(0,0,0,mi);      // pre-decay four-vector
  //  double masses[3] = {0,m_e,mass};  // neutrino, electron, nucleus
  //  rg.SetDecay(p0,3,masses);
  //  double wmax = rg.GetWtMax();
  //  for (int igen=0;igen<1000;igen++)   // cap the retries at 1000
  //    {
  //      double weight = rg.Generate();  // want to unweight these if possible
  //      TLorentzVector *e4v = rg.GetDecay(1);  // get electron four-vector
  //      p = e4v->P();
  //  if (weight >= wmax * flat.fire()) break;
  //   }
  //return p;
  //}




}//end namespace evgen

namespace evgen{

  DEFINE_ART_MODULE(RadioGen)

}//end namespace evgen

#endif
////////////////////////////////////////////////////////////////////////
