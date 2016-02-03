/**
 * @file   DumpMCParticles_module.cc
 * @brief  Module dumping MCarticle information on screen
 * @date   December 3rd, 2015
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 */


// nutools libraries
#include "SimulationBase/MCParticle.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


namespace sim {
  class DumpMCParticles;
} // namespace sim

namespace {
  using namespace fhicl;
  
  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> InputParticles {
      Name("InputParticles"),
      Comment("data product with the MC particles to be dumped")
      };
    
    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the output stream (managed by the message facility)"),
      "DumpMCParticles" /* default value */
      };
    
    fhicl::Atom<unsigned int> PointsPerLine {
      Name("PointsPerLine"),
      Comment("trajectory points printed per line (default: 2; 0 = skip them)"),
      2 /* default value */
      };
    
  }; // struct Config
  
  
} // local namespace


class sim::DumpMCParticles: public art::EDAnalyzer {
    public:
  // type to enable module parameters description by art
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  /// Configuration-checking constructor
  explicit DumpMCParticles(Parameters const& config);
  
  // Plugins should not be copied or assigned.
  DumpMCParticles(DumpMCParticles const&) = delete;
  DumpMCParticles(DumpMCParticles &&) = delete;
  DumpMCParticles& operator = (DumpMCParticles const&) = delete;
  DumpMCParticles& operator = (DumpMCParticles &&) = delete;
  
  
  // Operates on the event
  void analyze(art::Event const& event) override;
  
  
  /**
   * @brief Dumps the content of the specified particle in the output stream
   * @tparam Stream the type of output stream
   * @param out the output stream
   * @param particle the particle to be dumped
   * @param indent base indentation string (default: none)
   * @param bIndentFirst if first output line should be indented (default: yes)
   * 
   * The indent string is prepended to every line of output, with the possible
   * exception of the first one, in case bIndentFirst is true.
   * 
   * The output starts on the current line, and the last line is NOT broken.
   */
  template <typename Stream>
  void DumpMCParticle(
    Stream&& out, simb::MCParticle const& particle,
    std::string indent = "", bool bIndentFirst = true
    ) const;
  
  
    private:
  
  art::InputTag fInputParticles; ///< name of MCParticle's data product
  std::string fOutputCategory; ///< name of the stream for output
  unsigned int fPointsPerLine; ///< trajectory points per output line
  
}; // class sim::DumpMCParticles


//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
sim::DumpMCParticles::DumpMCParticles(Parameters const& config)
  : EDAnalyzer(config)
  , fInputParticles(config().InputParticles())
  , fOutputCategory(config().OutputCategory())
  , fPointsPerLine(config().PointsPerLine())
{}

//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpMCParticles::DumpMCParticle(
  Stream&& out, simb::MCParticle const& particle,
  std::string indent /* = "" */, bool bIndentFirst /* = true */
) const {
  if (bIndentFirst) out << indent;
  out << "ID=" << particle.TrackId() << "  PDG ID " << particle.PdgCode()
    << " mass=" << particle.Mass() << " GeV/c2 "
    << " status=" << particle.StatusCode();
  if (particle.Weight() != 1.0) out << " weight=" << particle.Weight();
  if (particle.Rescatter()) out << " rescattered";
  out
    << " at vertex (" << particle.Gvx() << ", " << particle.Gvy()
    << ", " << particle.Gvz() << "; " << particle.Gvt() << ")";
  const unsigned int nPoints = particle.NumberTrajectoryPoints();
  if (nPoints == 0) {
    out << " (no points!)";
  }
  else {
    TLorentzVector const& start = particle.Position();
    TLorentzVector const& start_mom = particle.Momentum();
    out << "\n" << indent << "created at ("
      << start.X() << ", " << start.Y() << ", " << start.Z() << "; "
      << start.T()
      << ") cm by ";
    if (particle.Mother() == 0) out << "the gods";
    else {
      out << (particle.Process().empty()? "magics": particle.Process())
        << " from ID=" << particle.Mother();
    }
    out << " with momentum ("
      << start_mom.X() << ", " << start_mom.Y() << ", " << start_mom.Z() << "; "
      << start_mom.T() << ") GeV/c";
    out << "\n" << indent;
    if (nPoints == 1) {
      out << "still alive!";
    }
    else {
      TLorentzVector const& stop = particle.EndPosition();
      TLorentzVector const& stop_mom = particle.EndMomentum();
      out << "stops at ("
        << stop.X() << ", " << stop.Y() << ", " << stop.Z() << "; "
        << stop.T() << ") cm with momentum ("
        << start_mom.X() << ", " << start_mom.Y()
        << ", " << start_mom.Z() << "; " << start_mom.T() << ") GeV/c by "
        << (particle.EndProcess().empty()? "magics": particle.EndProcess())
        << " into ";
      const unsigned int nDaughters = particle.NumberDaughters();
      if (nDaughters == 0) out << "nothing";
      else if (nDaughters == 1)
        out << "particle ID=" << particle.FirstDaughter();
      else {
        out << nDaughters << " particles from ID=" << particle.FirstDaughter()
          << " to ID=" << particle.LastDaughter();
      }
    } // if more than one point
    simb::MCTrajectory const& traj = particle.Trajectory();
    
    out << "\n" << indent << "comes with a trajectory " << traj.TotalLength()
      << " cm long in " << nPoints << " points";
    if (fPointsPerLine > 0) {
      out << ":";
      unsigned int page = 0;
      for (auto const& pair: traj) {
        if (page-- == 0) {
          out << "\n" << indent << "  ";
          page = fPointsPerLine - 1;
        }
        else out << " -- ";
        
        TLorentzVector const& pos = pair.first;
        out << "(" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << "; "
          << pos.T() << ")";
      } // for trajectory points
    } // if printing points
  } // if has points
  
} // sim::DumpMCParticles::DumpMCParticle()


//------------------------------------------------------------------------------
void sim::DumpMCParticles::analyze(art::Event const& event) {
  
  // get the particles from the event
  auto const& Particles
    = *(event.getValidHandle<std::vector<simb::MCParticle>>(fInputParticles));
  
  mf::LogInfo(fOutputCategory)
    << "Data product '" << fInputParticles.encode() << "' contains "
    << Particles.size() << " MCParticle's";
  
  unsigned int iParticle = 0;
  mf::LogVerbatim log(fOutputCategory);
  for (simb::MCParticle const& particle: Particles) {
    
    // a bit of a header
    log << "\n[#" << (iParticle++) << "] ";
    DumpMCParticle(log, particle, "  ", false);
    
  } // for
  
} // sim::DumpMCParticles::analyze()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpMCParticles)

//------------------------------------------------------------------------------
