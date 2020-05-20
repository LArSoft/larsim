////////////////////////////////////////////////////////////////////////
// Class:       NeutronOsc
// Module Type: producer
// GENIE neutron-antineutron oscillation generator
//
// Adapted from NucleonDecay_module.cc (tjyang@fnal.gov)
// by jhewes15@fnal.gov
//
//  Neutron-antineutron oscillation mode ID:
// ---------------------------------------------------------
//  ID |   Decay Mode
//     |
// ---------------------------------------------------------
//   0 |    Random oscillation mode
//   1 |    p + nbar --> \pi^{+} + \pi^{0}
//   2 |    p + nbar --> \pi^{+} + 2\pi^{0}
//   3 |    p + nbar --> \pi^{+} + 3\pi^{0}
//   4 |    p + nbar --> 2\pi^{+} + \pi^{-} + \pi^{0}
//   5 |    p + nbar --> 2\pi^{+} + \pi^{-} + 2\pi^{0}
//   6 |    p + nbar --> 2\pi^{+} + \pi^{-} + 2\omega^{0}
//   7 |    p + nbar --> 3\pi^{+} + 2\pi^{-} + \pi^{0}
//   8 |    n + nbar --> \pi^{+} + \pi^{-}
//   9 |    n + nbar --> 2\pi^{0}
//  10 |    n + nbar --> \pi^{+} + \pi^{-} + \pi^{0}
//  11 |    n + nbar --> \pi^{+} + \pi^{-} + 2\pi^{0}
//  12 |    n + nbar --> \pi^{+} + \pi^{-} + 3\pi^{0}
//  13 |    n + nbar --> 2\pi^{+} + 2\pi^{-}
//  14 |    n + nbar --> 2\pi^{+} + 2\pi^{-} + \pi^{0}
//  15 |    n + nbar --> \pi^{+} + \pi^{-} + \omega^{0}
//  16 |    n + nbar --> 2\pi^{+} + 2\pi^{-} + 2\pi^{0}
// ---------------------------------------------------------
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// GENIE includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/EventGen/EventRecord.h"
#include "Physics/NNBarOscillation/NNBarOscMode.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Utils/AppInit.h"

// larsoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nurandom/RandomUtils/NuRandomService.h"

// c++ includes
#include <memory>
#include <string>

#include "CLHEP/Random/RandFlat.h"

namespace evgen {
  class NeutronOsc;
}

class evgen::NeutronOsc : public art::EDProducer {
public:
  explicit NeutronOsc(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutronOsc(NeutronOsc const &) = delete;
  NeutronOsc(NeutronOsc &&) = delete;
  NeutronOsc & operator = (NeutronOsc const &) = delete;
  NeutronOsc & operator = (NeutronOsc &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginRun(art::Run& run) override;

private:

  // Additional functions
  int SelectAnnihilationMode(int pdg_code);

  // Declare member data here.
  const genie::EventRecordVisitorI * mcgen;
  genie::NNBarOscMode_t gOptDecayMode    = genie::kNONull;             // neutron-antineutron oscillation mode
  CLHEP::RandFlat flatDist;

};


evgen::NeutronOsc::NeutronOsc(fhicl::ParameterSet const & p)
  : art::EDProducer{p}
  // create a default random engine; obtain the random seed from NuRandomService,
  // unless overridden in configuration with key "Seed"
  , flatDist{art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, p, "Seed")}
{
  genie::PDGLibrary::Instance(); //Ensure Messenger is started first in GENIE.

  string sname   = "genie::EventGenerator";
  // GENIE v2 // string sconfig = "NeutronOsc";
  string sconfig = "NNBarOsc";
  // GENIE v3 needs a tune (even if irrelevant)
  evgb::SetEventGeneratorListAndTune("Default","Default");

  genie::AlgFactory * algf = genie::AlgFactory::Instance();
  mcgen =
    dynamic_cast<const genie::EventRecordVisitorI *> (algf->GetAlgorithm(sname,sconfig));
  if(!mcgen) {
    throw cet::exception("NeutronOsc") << "Couldn't instantiate the neutron-antineutron oscillation generator";
  }
  int fDecayMode = p.get<int>("DecayMode");
  gOptDecayMode = (genie::NNBarOscMode_t) fDecayMode;

  produces< std::vector<simb::MCTruth> >();
  produces< sumdata::RunData, art::InRun >();

  unsigned int seed = art::ServiceHandle<rndm::NuRandomService>()->getSeed();
  genie::utils::app_init::RandGen(seed);
}

void evgen::NeutronOsc::produce(art::Event & e)
{
  // Implementation of required member function here.
  genie::EventRecord * event = new genie::EventRecord;
  int target = 1000180400;  //Only use argon target
  int decay  = SelectAnnihilationMode(target);
  genie::Interaction * interaction = genie::Interaction::NOsc(target,decay);
  event->AttachSummary(interaction);

  // Simulate decay
  mcgen->ProcessEventRecord(event);

//  genie::Interaction *inter = event->Summary();
//  const genie::InitialState &initState = inter->InitState();
//  std::cout<<"initState = "<<initState.AsString()<<std::endl;
//  const genie::ProcessInfo &procInfo = inter->ProcInfo();
//  std::cout<<"procInfo = "<<procInfo.AsString()<<std::endl;
  MF_LOG_DEBUG("NeutronOsc")
    << "Generated event: " << *event;

  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
  simb::MCTruth truth;

  art::ServiceHandle<geo::Geometry const> geo;

  // Find boundary of active volume
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (size_t i = 0; i<geo->NTPC(); ++i){
    const geo::TPCGeo &tpc = geo->TPC(i);
    if (minx>tpc.MinX()) minx = tpc.MinX();
    if (maxx<tpc.MaxX()) maxx = tpc.MaxX();
    if (miny>tpc.MinY()) miny = tpc.MinY();
    if (maxy<tpc.MaxY()) maxy = tpc.MaxY();
    if (minz>tpc.MinZ()) minz = tpc.MinZ();
    if (maxz<tpc.MaxZ()) maxz = tpc.MaxZ();
  }

  // Assign vertice position
  double X0 = flatDist.fire( minx, maxx );
  double Y0 = flatDist.fire( miny, maxy );
  double Z0 = flatDist.fire( minz, maxz );

  TIter partitr(event);
  genie::GHepParticle *part = 0;
  // GHepParticles return units of GeV/c for p.  the V_i are all in fermis
  // and are relative to the center of the struck nucleus.
  // add the vertex X/Y/Z to the V_i for status codes 0 and 1
  int trackid = 0;
  std::string primary("primary");

  while( (part = dynamic_cast<genie::GHepParticle *>(partitr.Next())) ){

    simb::MCParticle tpart(trackid,
                           part->Pdg(),
                           primary,
                           part->FirstMother(),
                           part->Mass(),
                           part->Status());

    TLorentzVector pos(X0, Y0, Z0, 0);
    TLorentzVector mom(part->Px(), part->Py(), part->Pz(), part->E());
    tpart.AddTrajectoryPoint(pos,mom);
    if(part->PolzIsSet()) {
      TVector3 polz;
      part->GetPolarization(polz);
      tpart.SetPolarization(polz);
    }
    truth.Add(tpart);

    ++trackid;
  }// end loop to convert GHepParticles to MCParticles
  truth.SetOrigin(simb::kUnknown);
  truthcol->push_back(truth);
  //FillHistograms(truth);
  e.put(std::move(truthcol));

  delete event;
}

void evgen::NeutronOsc::beginRun(art::Run& run)
{
  art::ServiceHandle<geo::Geometry const> geo;
  run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
}


int evgen::NeutronOsc::SelectAnnihilationMode(int pdg_code)
{
  // if the mode is set to 'random' (the default), pick one at random!
  if ((int)gOptDecayMode == 0) {
    int mode;

    std::string pdg_string = std::to_string(static_cast<long long>(pdg_code));
    if (pdg_string.size() != 10) {
      std::cout << "Expecting PDG code to be a 10-digit integer; instead, it's the following: " << pdg_string << std::endl;
      exit(1);
    }

    // count number of protons & neutrons
    int n_nucleons = std::stoi(pdg_string.substr(6,3)) - 1;
    int n_protons  = std::stoi(pdg_string.substr(3,3));

    // factor proton / neutron ratio into branching ratios
    double proton_frac  = ((double)n_protons) / ((double)n_nucleons);
    double neutron_frac = 1 - proton_frac;

    // set branching ratios, taken from bubble chamber data
    const int n_modes = 16;
    double br [n_modes] = { 0.010, 0.080, 0.100, 0.220,
                            0.360, 0.160, 0.070, 0.020,
                            0.015, 0.065, 0.110, 0.280,
                            0.070, 0.240, 0.100, 0.100 };

    for (int i = 0; i < n_modes; i++) {
      if (i < 7)
        br[i] *= proton_frac;
      else
        br[i] *= neutron_frac;
    }

    // randomly generate a number between 1 and 0
    double p = flatDist.fire();

    // loop through all modes, figure out which one our random number corresponds to
    double threshold = 0;
    for (int i = 0; i < n_modes; i++) {
      threshold += br[i];
      if (p < threshold) {
        // once we've found our mode, return it!
        mode = i + 1;
        return mode;
      }
    }

    // error message, in case the random number selection fails
    std::cout << "Random selection of final state failed!" << std::endl;
    exit(1);
  }

  // if specific annihilation mode specified, just use that
  else {
    int mode = (int) gOptDecayMode;
    return mode;
  }
}

DEFINE_ART_MODULE(evgen::NeutronOsc)