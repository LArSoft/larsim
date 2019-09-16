/**
 * @file   DumpMCTracks_module.cc
 * @brief  Module dumping MCTrack information on screen
 * @date   March 21st, 2016
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 */

// LArSoft libraries
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/MCBase/MCTrack.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <string>
#include <algorithm> // std::max()


namespace sim {
  class DumpMCTracks;
} // namespace sim

namespace {
  using namespace fhicl;

  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> InputTracks {
      Name("InputTracks"),
      Comment("data product with the MC tracks to be dumped")
      };

    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the output stream (managed by the message facility)"),
      "DumpMCTracks" /* default value */
      };

  }; // struct Config


  /// Returns a string describing the type of origin
  std::string OriginDescription(simb::Origin_t origin) {
    switch (origin) {
      case simb::kUnknown:           return "unknown";
      case simb::kBeamNeutrino:      return "beam neutrino";
      case simb::kCosmicRay:         return "cosmic rays";
      case simb::kSuperNovaNeutrino: return "supernova neutrinos";
      case simb::kSingleParticle:    return "single particles thrown at the detector";
      // default case is deliberately missing,
      // so that on new additions in nutools the compiler will complain
    } // switch
    throw art::Exception(art::errors::LogicError)
      << "Unexpected origin type #" << ((int) origin) << "\n";
  } // OriginDescription()


  /// Prints a MC step into a stream
  template <typename Stream>
  void PrintMCStep(Stream&& out, sim::MCStep const& step) {
    out << "("
      << step.X() << ", " << step.Y() << ", " << step.Z() << ") cm, t="
      << step.T() << " ns; momentum ("
      << step.Px() << ", " << step.Py() << ", " << step.Pz() << "; "
      << step.E() << ") MeV/c";
  } // PrintMCStep()


} // local namespace


class sim::DumpMCTracks: public art::EDAnalyzer {
    public:
  // type to enable module parameters description by art
  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Configuration-checking constructor
  explicit DumpMCTracks(Parameters const& config);

  // Plugins should not be copied or assigned.
  DumpMCTracks(DumpMCTracks const&) = delete;
  DumpMCTracks(DumpMCTracks &&) = delete;
  DumpMCTracks& operator = (DumpMCTracks const&) = delete;
  DumpMCTracks& operator = (DumpMCTracks &&) = delete;


  // Operates on the event
  virtual void analyze(art::Event const& event) override;


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
  void DumpMCTrack(
    Stream&& out, sim::MCTrack const& track,
    std::string indent = "", bool bIndentFirst = true
    ) const;


    private:

  art::InputTag fInputTracks; ///< name of MCTrack's data product
  std::string fOutputCategory; ///< name of the stream for output

}; // class sim::DumpMCTracks


//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
sim::DumpMCTracks::DumpMCTracks(Parameters const& config)
  : EDAnalyzer(config)
  , fInputTracks(config().InputTracks())
  , fOutputCategory(config().OutputCategory())
{}

//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpMCTracks::DumpMCTrack(
  Stream&& out, sim::MCTrack const& track,
  std::string indent /* = "" */, bool bIndentFirst /* = true */
) const {
  if (bIndentFirst) out << indent;
  out
    << "from GEANT track ID=" << track.TrackID()
    << " PDG ID=" << track.PdgCode()
    << " from " << OriginDescription(track.Origin())
    << " via '" << track.Process() << "'";
  out << "\n" << indent
    << "  starting at ";
  ::PrintMCStep(out, track.Start());
  out << "\n" << indent
    << "  ending at ";
  ::PrintMCStep(out, track.End());

  std::vector<std::vector<double>> const& dQdx = track.dQdx(); // dQdx[MCStep][plane]
  std::vector<double> const& dEdx = track.dEdx(); // dEdx[MCStep]
  size_t const nQSteps = dQdx.size(), nESteps = dEdx.size();
  size_t const nSteps = std::max(nQSteps, nESteps);
  out << "\n" << indent;
  if (nSteps > 0) {
    out
      << "energy information for " << nSteps
      << " steps (dE/dX in MeV/cm, then dQ/dx per plane):";
    for (size_t iStep = 0; iStep < nSteps; ++iStep) {
      out << "\n" << indent
        << "  [#" << iStep << "] dE/dx=";
      if (iStep < nESteps) out << dEdx[iStep];
      else                 out << "<N/A>";
      out << "; dQ/dx:";
      if (iStep < nQSteps) {
        std::vector<double> const& step_dQdx = dQdx[iStep]; // dQdx[plane]
        for (size_t iPlane = 0; iPlane < step_dQdx.size(); ++iPlane) {
          out << " [#" << iPlane << "] " << step_dQdx[iPlane];
        } // for plane
      }
      else                 out << "<N/A>";
    } // for iStep
  }
  else out << "no energy or charge information available";

  out << "\n" << indent
    << "mother ID=" << track.MotherTrackID()
    << " PDG ID=" << track.MotherPdgCode()
    << " via '" << track.MotherProcess() << "'";
  out << "\n" << indent
    << "  starting at ";
  ::PrintMCStep(out, track.MotherStart());
  out << "\n" << indent
    << "  ending at ";
  ::PrintMCStep(out, track.MotherEnd());

  out << "\n" << indent
    << "ancestor ID=" << track.AncestorTrackID()
    << " PDG ID=" << track.AncestorPdgCode()
    << " via '" << track.AncestorProcess() << "'";
  out << "\n" << indent
    << "  starting at ";
  ::PrintMCStep(out, track.AncestorStart());
  out << "\n" << indent
    << "  ending at ";
  ::PrintMCStep(out, track.AncestorEnd());

} // sim::DumpMCTracks::DumpMCTrack()


//------------------------------------------------------------------------------
void sim::DumpMCTracks::analyze(art::Event const& event) {

  // get the particles from the event
  auto const& Tracks
    = *(event.getValidHandle<std::vector<sim::MCTrack>>(fInputTracks));

  mf::LogVerbatim(fOutputCategory)
    << "Event " << event.id() << ": data product '"
    << fInputTracks.encode() << "' contains "
    << Tracks.size() << " MCTrack objects";

  unsigned int iTrack = 0;
  mf::LogVerbatim log(fOutputCategory);
  for (sim::MCTrack const& track: Tracks) {

    // a bit of a header
    log << "\n[#" << (iTrack++) << "] ";
    DumpMCTrack(log, track, "  ", false);

  } // for
  log << "\n";

} // sim::DumpMCTracks::analyze()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpMCTracks)

//------------------------------------------------------------------------------
