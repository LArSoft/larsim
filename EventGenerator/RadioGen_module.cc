////////////////////////////////////////////////////////////////////////
/// \file  RadioGen_module.cc
/// \brief Generator for radiological decays
///
/// Module designed to produce a set list of particles for MC to model radiological decays
///
/// \version $Id: RadioGen_module.cc,v 1.0 2014/09/05 15:02:00 trj Exp $
/// \author  trj@fnal.gov
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

// nutools includes
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "EventGeneratorBase/evgenbase.h"

// lar includes
#include "Geometry/Geometry.h"
#include "SummaryData/RunData.h"

// root includes

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGenPhaseSpace.h"

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

    double betaphasespace(double mass, double q);

    unsigned int        fSeed;           ///< random number seed    
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

    const double gevperamu = 0.931494061;
    const double m_e = 0.000511;  // mass of electron in GeV

    TGenPhaseSpace rg;  // put this here so we don't constantly construct and destruct it

  };
}

namespace evgen{

  //____________________________________________________________________________
  RadioGen::RadioGen(fhicl::ParameterSet const& pset)
  {

    this->reconfigure(pset);

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    fSeed = pset.get< unsigned int >("Seed", evgb::GetRandomNumberSeed());

    createEngine( fSeed );

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
    // do not put fSeed in reconfigure because we don't want to reset 
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

	int pdgid=0;  // electron=11, photon=22
	double p = 0; // generated momentum (GeV)
	double m = 0; // mass of daughter particle

	if (fNuclide[i] == "39Ar")
	  {
	    double mass = 38.964313 * 0.931494061 - 18*m_e; // 39Ar nucleus mass in GeV.
	    double q = 0.000565;
	    p = betaphasespace(mass,q); // first attempt -- use phase space
	    pdgid = 11;  // generate beta spectrum
	    m = m_e;
	  }
	else
	  {
	    throw cet::exception("RadioGen") << "Unimplemented nuclide: " << fNuclide[i] << "\n";
	  }

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
	//pvec.Print();
 
	// set track id to a negative serial number as these are all primary particles and have id <= 0
	int trackid = trackidcounter;
	trackidcounter--;
	std::string primary("primary");
	simb::MCParticle part(trackid, pdgid, primary);
	part.AddTrajectoryPoint(pos, pvec);
	mct.Add(part);
      }
  }

  // phase space generator for beta decay

  double RadioGen::betaphasespace(double mass, double q)
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat     flat(engine);
    double p = 0;
    double mi = mass+q+m_e;
    TLorentzVector p0(0,0,0,mi);      // pre-decay four-vector
    double masses[3] = {0,m_e,mass};  // neutrino, electron, nucleus
    rg.SetDecay(p0,3,masses);
    double wmax = rg.GetWtMax();
    for (int igen=0;igen<1000;igen++)   // cap the retries at 1000
      {
        double weight = rg.Generate();  // want to unweight these if possible
        TLorentzVector *e4v = rg.GetDecay(1);  // get electron four-vector
        p = e4v->P();
	if (weight >= wmax * flat.fire()) break;
      }
    return p;
  }

}//end namespace evgen

namespace evgen{

  DEFINE_ART_MODULE(RadioGen)

}//end namespace evgen

#endif
////////////////////////////////////////////////////////////////////////
