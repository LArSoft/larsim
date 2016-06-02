////////////////////////////////////////////////////////////////////////
/// \file  RadioGen_module.cc
/// \brief Generator for radiological decays
///
/// Module designed to produce a set list of particles for MC to model radiological decays
///
/// \version $Id: RadioGen_module.cc,v 1.0 2014/09/05 15:02:00 trj Exp $
/// \author  trj@fnal.gov
//           Rn222 generation feature added by gleb.sinev@duke.edu 
//           (based on a generator by jason.stock@mines.sdsmt.edu)
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_RADIOLOGICAL
#define EVGEN_RADIOLOGICAL

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <memory>
#include <iterator>
#include <vector>
#include <sys/stat.h>

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
#include "cetlib/search_path.h"

// art extensions
#include "larsim/RandomUtils/LArSeedService.h"

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

    void SampleOne(unsigned int   i, 
		   simb::MCTruth &mct);        

    void readfile(std::string nuclide, std::string filename);
    void samplespectrum(std::string nuclide, int &itype, double &t, double &m, double &p);  

    // recoded so as to use the LArSoft-managed random number generator
    double samplefromth1d(TH1D *hist);

    // itype = pdg code:  1000020040: alpha.  itype=11: beta. -11: positron,  itype=22: gamma.  -1: error
    // t is the kinetic energy in GeV  (= E-m)
    // m = mass in GeV
    // p = momentum in GeV

    //double betaphasespace(double mass, double q); // older parameterization.

    std::vector<std::string> fNuclide;   ///< List of nuclides to simulate.  Example:  "39Ar".
    std::vector<double> fBq;             ///< Radioactivity in Becquerels (decay per sec) per cubic cm.
    std::vector<double> fT0;             ///< Beginning of time window to simulate in ns
    std::vector<double> fT1;             ///< End of time window to simulate in ns
    std::vector<double> fX0;             ///< Bottom corner x position (cm) in world coordinates 
    std::vector<double> fY0;             ///< Bottom corner y position (cm) in world coordinates
    std::vector<double> fZ0;             ///< Bottom corner z position (cm) in world coordinates
    std::vector<double> fX1;             ///< Top corner x position (cm) in world coordinates 
    std::vector<double> fY1;             ///< Top corner y position (cm) in world coordinates
    std::vector<double> fZ1;             ///< Top corner z position (cm) in world coordinates
    int trackidcounter;                  ///< Serial number for the MC track ID

    // leftovers from the phase space generator
    // const double gevperamu = 0.931494061;
    // TGenPhaseSpace rg;  // put this here so we don't constantly construct and destruct it

    const double m_e = 0.000510998928;  // mass of electron in GeV
    const double m_alpha = 3.727379240; // mass of an alpha particle in GeV

    std::vector<std::string> spectrumname;
    std::vector<TH1D*> alphaspectrum;
    std::vector<double> alphaintegral;
    std::vector<TH1D*> betaspectrum;
    std::vector<double> betaintegral;
    std::vector<TH1D*> gammaspectrum;
    std::vector<double> gammaintegral;

  };
}

namespace evgen{

  //____________________________________________________________________________
  RadioGen::RadioGen(fhicl::ParameterSet const& pset)
  {

    this->reconfigure(pset);

    // create a default random engine; obtain the random seed from LArSeedService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<sim::LArSeedService>()
      ->createEngine(*this, pset, "Seed");

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
    fBq            = p.get< std::vector<double> >("BqPercc");
    fT0            = p.get< std::vector<double> >("T0");
    fT1            = p.get< std::vector<double> >("T1");
    fX0            = p.get< std::vector<double> >("X0");
    fY0            = p.get< std::vector<double> >("Y0");
    fZ0            = p.get< std::vector<double> >("Z0");
    fX1            = p.get< std::vector<double> >("X1");
    fY1            = p.get< std::vector<double> >("Y1");
    fZ1            = p.get< std::vector<double> >("Z1");

    // check for consistency of vector sizes
    
    unsigned int nsize = fNuclide.size();
    if (  fBq.size() != nsize ) throw cet::exception("RadioGen") << "Different size Bq vector and Nuclide vector\n";
    if (  fT0.size() != nsize ) throw cet::exception("RadioGen") << "Different size T0 vector and Nuclide vector\n";
    if (  fT1.size() != nsize ) throw cet::exception("RadioGen") << "Different size T1 vector and Nuclide vector\n";
    if (  fX0.size() != nsize ) throw cet::exception("RadioGen") << "Different size X0 vector and Nuclide vector\n";
    if (  fY0.size() != nsize ) throw cet::exception("RadioGen") << "Different size Y0 vector and Nuclide vector\n";
    if (  fZ0.size() != nsize ) throw cet::exception("RadioGen") << "Different size Z0 vector and Nuclide vector\n";
    if (  fX1.size() != nsize ) throw cet::exception("RadioGen") << "Different size X1 vector and Nuclide vector\n";
    if (  fY1.size() != nsize ) throw cet::exception("RadioGen") << "Different size Y1 vector and Nuclide vector\n";
    if (  fZ1.size() != nsize ) throw cet::exception("RadioGen") << "Different size Z1 vector and Nuclide vector\n";

    readfile("39Ar","Argon_39.root");
    readfile("60Co","Cobalt_60.root");
    readfile("85Kr","Krypton_85.root");
    readfile("40K","Potassium_40.root");
    readfile("232Th","Thorium_232.root");
    readfile("238U","Uranium_238.root");

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

    LOG_DEBUG("RadioGen") << truth;
    truthcol->push_back(truth);
    evt.put(std::move(truthcol));
    return;
  }

  //____________________________________________________________________________
  // Generate radioactive decays per nuclide per volume according to the FCL parameters

  void RadioGen::SampleOne(unsigned int i, simb::MCTruth &mct){

    // get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat     flat(engine);
    CLHEP::RandPoisson  poisson(engine);

    // figure out how many decays to generate

    double rate = fabs( fBq[i] * (fT1[i] - fT0[i]) * (fX1[i] - fX0[i]) * (fY1[i] - fY0[i]) * (fZ1[i] - fZ0[i]) ) / 1.0E9;
    long ndecays = poisson.shoot(rate);

    for (unsigned int idecay=0; idecay<ndecays; idecay++)
      {
	// generate just one particle at a time.  Need a little recoding if a radioactive
	// decay generates multiple daughters that need simulation

	int pdgid=0;  // electron=11, photon=22, alpha = 1000020040
        double t = 0; // kinetic energy of particle GeV
	double m = 0; // mass of daughter particle GeV
	double p = 0; // generated momentum (GeV)

        // Treat 222Rn separately 
        if (fNuclide[i] == "222Rn")
        {
          pdgid = 1000020040;
          t     = 0.00548952;
          m     = m_alpha;
          double energy = t + m;
          double p2     = energy*energy - m*m;
          if (p2 > 0) p = TMath::Sqrt(p2);
          else        p = 0;
        }
        else samplespectrum(fNuclide[i],pdgid,t,m,p);


	// uniformly distributed in position and time

	TLorentzVector pos( fX0[i] + flat.fire()*(fX1[i] - fX0[i]),
			    fY0[i] + flat.fire()*(fY1[i] - fY0[i]),
			    fZ0[i] + flat.fire()*(fZ1[i] - fZ0[i]),
			    fT0[i] + flat.fire()*(fT1[i] - fT0[i]) );

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

	// set track id to a negative serial number as these are all primary particles and have id <= 0
	int trackid = trackidcounter;
	trackidcounter--;
	std::string primary("primary");
	// alpha particles need a little help since they're not in the TDatabasePDG table
	// so don't rely so heavily on default arguments to the MCParticle constructor
	if (pdgid == 1000020040)
	  {
	    simb::MCParticle part(trackid, pdgid, primary,-1,m,1);
	    part.AddTrajectoryPoint(pos, pvec);
	    mct.Add(part);
	  }
	else
	  {
	    simb::MCParticle part(trackid, pdgid, primary);
	    part.AddTrajectoryPoint(pos, pvec);
	    mct.Add(part);
	  }
      }
  }

  // only reads those files that are on the fNuclide list.  Copy information from the TGraphs to TH1D's

  void RadioGen::readfile(std::string nuclide, std::string filename)
  {
    int ifound = 0;
    for (size_t i=0; i<fNuclide.size(); i++)
      {
	if (fNuclide[i] == nuclide)
	  {
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

    if (alphagraph)
      {
	int np = alphagraph->GetN();
	double *y = alphagraph->GetY();
	std::string name;
	name = "RadioGen_";
	name += nuclide;
	name += "_Alpha";
	TH1D *alphahist = (TH1D*) new TH1D(name.c_str(),"Alpha Spectrum",np,0,np-1);
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
	TH1D *betahist = (TH1D*) new TH1D(name.c_str(),"Beta Spectrum",np,0,np-1);

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
	TH1D *gammahist = (TH1D*) new TH1D(name.c_str(),"Gamma Spectrum",np,0,np-1);
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

    f.Close();
    TH1::AddDirectory(addStatus);

    double total = alphaintegral.back() + betaintegral.back() + gammaintegral.back();
    if (total>0)
      {
	alphaintegral.back() /= total;
	betaintegral.back() /= total;
	gammaintegral.back() /= total;
      }
  }


  void RadioGen::samplespectrum(std::string nuclide, int &itype, double &t, double &m, double &p)
  {

    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat  flat(engine);

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
	else if ( gammaspectrum[inuc] != 0)
	  {
	    itype = 22; // gamma
	    m = 0;
	    t = samplefromth1d(gammaspectrum[inuc])/1000000.0;
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
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat  flat(engine);

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
  //	if (weight >= wmax * flat.fire()) break;
  //   }
  //return p;
  //}




}//end namespace evgen

namespace evgen{

  DEFINE_ART_MODULE(RadioGen)

}//end namespace evgen

#endif
////////////////////////////////////////////////////////////////////////
