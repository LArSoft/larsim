/**
 * @file TextFileGen_module.cc
 * @brief Producer generating Monte Carlo truth record in LArSoft format from a text file
 * @date Mon Apr  8 09:20:02 2013
 * @author Brian Rebel
 */
/**
 * @class evgen::TextFileGen
 *
 *  This module assumes that the input file has the hepevt format for
 *  each event to be simulated.  See
 *
 *  http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node39.html
 *
 *  for details on the format.  In brief each event contains at least two
 *  lines.  The first line contains two entries, the event number (which is
 *  ignored in ART/LArSoft) and the number of particles in the event.  Each
 *  following line containes 15 entries to describe each particle.  The entries
 *  are:
 *
 *  1.  status code (should be set to 1 for any particle to be tracked, others
 *      won't be tracked)
 *  2.  the pdg code for the particle
 *  3.  the entry of the first mother for this particle in the event,
 *      0 means no mother
 *  4.  the entry of the second mother for this particle in the event,
 *      0 means no mother
 *  5. the entry of the first daughter for this particle in the event,
 *      0 means no daughter
 *  6. the entry of the second daughter for this particle in the event,
 *      0 means no daughter
 *  7. x component of the particle momentum
 *  8. y component of the particle momentum
 *  9. z component of the particle momentum
 *  10. energy of the particle
 *  11. mass of the particle
 *  12. x position of the particle initial position
 *  13. y position of the particle initial position
 *  14. z position of the particle initial position
 *  15. time of the particle production
 *
 *  For example, if you want to simulate a single muon with a 5 GeV energy
 *  moving only in the z direction, the entry would be
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  0 1
 *  1 13 0 0 0 0 0. 0. 1.0 5.0011 0.105 1.0 1.0 1.0 0.0
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  There are some assumptions that go into using this format that may not
 *  be obvious.  The first is that only particles with status code = 1
 *  are tracked in the LArSoft/Geant4 combination making the mother daughter
 *  relations somewhat irrelevant.  That also means that you should let
 *  Geant4 handle any decays.
 *
 *  The units in LArSoft are cm for distances and ns for time.
 *  The use of `TLorentzVector` below does not imply space and time have the same units
 *   (do not use `TLorentzVector::Boost()`).
 */
#include <fstream>
#include <memory>
#include <string>
#include <utility>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TLorentzVector.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace evgen {
  class TextFileGen;
}

class evgen::TextFileGen : public art::EDProducer {
public:
  explicit TextFileGen(fhicl::ParameterSet const& p);

  void produce(art::Event& e) override;
  void beginJob() override;
  void beginRun(art::Run& run) override;

private:
  std::pair<unsigned, unsigned> readEventInfo(std::istream& is);
  simb::MCTruth readNextHepEvt();
  unsigned long int fOffset;
  std::ifstream* fInputFile;
  std::string fInputFileName; ///< Name of text file containing events to simulate
  double fMoveY;              ///< Project particles to a new y plane.
};

//------------------------------------------------------------------------------
evgen::TextFileGen::TextFileGen(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fOffset{p.get<unsigned long int>("Offset")}
  , fInputFile(0)
  , fInputFileName{p.get<std::string>("InputFileName")}
  , fMoveY{p.get<double>("MoveY", -1e9)}
{
  if (fMoveY > -1e8) {
    mf::LogWarning("TextFileGen") << "Particles will be moved to a new plane y = " << fMoveY
                                  << " cm.\n";
  }

  produces<std::vector<simb::MCTruth>>();
  produces<sumdata::RunData, art::InRun>();
}

//------------------------------------------------------------------------------
void evgen::TextFileGen::beginJob()
{
  fInputFile = new std::ifstream(fInputFileName.c_str());

  // check that the file is a good one
  if (!fInputFile->good())
    throw cet::exception("TextFileGen")
      << "input text file " << fInputFileName << " cannot be read.\n";

  for (unsigned i = 0; i != fOffset; ++i) {
    auto const [eventNo, nparticles] = readEventInfo(*fInputFile);
    for (unsigned p = 0; p != nparticles; ++p) {
      constexpr auto all_chars_until = std::numeric_limits<unsigned>::max();
      fInputFile->ignore(all_chars_until, '\n');
    }
  }
}

//------------------------------------------------------------------------------
void evgen::TextFileGen::beginRun(art::Run& run)
{
  art::ServiceHandle<geo::Geometry const> geo;
  run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
}

//------------------------------------------------------------------------------
void evgen::TextFileGen::produce(art::Event& e)
{
  // check that the file is still good
  if (!fInputFile->good())
    throw cet::exception("TextFileGen")
      << "input text file " << fInputFileName << " cannot be read in produce().\n";

  //Now, read the Event to be used.

  // check that the file is still good
  if (!fInputFile->good())
    throw cet::exception("TextFileGen")
      << "input text file " << fInputFileName << " cannot be read in produce().\n";
  auto truthcol = std::make_unique<std::vector<simb::MCTruth>>();
  truthcol->push_back(readNextHepEvt());

  e.put(std::move(truthcol));
}

simb::MCTruth evgen::TextFileGen::readNextHepEvt()
{

  // declare the variables for reading in the event record
  int status = 0;
  int pdg = 0;
  int firstMother = 0;
  int secondMother = 0;
  int firstDaughter = 0;
  int secondDaughter = 0;
  double xMomentum = 0.;
  double yMomentum = 0.;
  double zMomentum = 0.;
  double energy = 0.;
  double mass = 0.;
  double xPosition = 0.;
  double yPosition = 0.;
  double zPosition = 0.;
  double time = 0.;

  // read in line to get event number and number of particles
  std::string oneLine;
  std::istringstream inputLine;
  simb::MCTruth nextEvent;
  auto const [eventNo, nParticles] = readEventInfo(*fInputFile);

  // now read in all the lines for the particles
  // in this interaction. only particles with
  // status = 1 get tracked in Geant4.
  for (unsigned short i = 0; i < nParticles; ++i) {
    std::getline(*fInputFile, oneLine);
    inputLine.clear();
    inputLine.str(oneLine);

    inputLine >> status >> pdg >> firstMother >> secondMother >> firstDaughter >> secondDaughter >>
      xMomentum >> yMomentum >> zMomentum >> energy >> mass >> xPosition >> yPosition >>
      zPosition >> time;

    //Project the particle to a new y plane
    if (fMoveY > -1e8) {
      double totmom = sqrt(pow(xMomentum, 2) + pow(yMomentum, 2) + pow(zMomentum, 2));
      double kx = xMomentum / totmom;
      double ky = yMomentum / totmom;
      double kz = zMomentum / totmom;
      if (ky) {
        double l = (fMoveY - yPosition) / ky;
        xPosition += kx * l;
        yPosition += ky * l;
        zPosition += kz * l;
      }
    }

    TLorentzVector pos(xPosition, yPosition, zPosition, time);
    TLorentzVector mom(xMomentum, yMomentum, zMomentum, energy);

    simb::MCParticle part(i, pdg, "primary", firstMother, mass, status);
    part.AddTrajectoryPoint(pos, mom);

    nextEvent.Add(part);

  } //  end loop on particles.

  return nextEvent;
}

std::pair<unsigned, unsigned> evgen::TextFileGen::readEventInfo(std::istream& iss)
{
  std::string line;
  getline(iss, line);
  std::istringstream buffer{line};

  // Parse read line for the event number and particles per event
  unsigned event, nparticles;
  buffer >> event >> nparticles;
  return {event, nparticles};
}

DEFINE_ART_MODULE(evgen::TextFileGen)
