/**
 * @file   DumpMCShowers_module.cc
 * @brief  Module dumping MCShower information on screen
 * @date   March 21st, 2016
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 */

// LArSoft libraries
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/MCBase/MCShower.h"

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
#include <iomanip>
#include <algorithm> // std::max()


namespace sim {
  class DumpMCShowers;
} // namespace sim

namespace {
  using namespace fhicl;

  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> InputShowers {
      Name("InputShowers"),
      Comment("data product with the MC showers to be dumped")
      };

    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the output stream (managed by the message facility)"),
      "DumpMCShowers" /* default value */
      };

    fhicl::Atom<unsigned int> DaughtersPerLine {
      Name("DaughtersPerLine"),
      Comment("daughter IDs to print on each output line"),
      12 /* default value */
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


class sim::DumpMCShowers: public art::EDAnalyzer {
    public:
  // type to enable module parameters description by art
  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Configuration-checking constructor
  explicit DumpMCShowers(Parameters const& config);

  // Plugins should not be copied or assigned.
  DumpMCShowers(DumpMCShowers const&) = delete;
  DumpMCShowers(DumpMCShowers &&) = delete;
  DumpMCShowers& operator = (DumpMCShowers const&) = delete;
  DumpMCShowers& operator = (DumpMCShowers &&) = delete;


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
  void DumpMCShower(
    Stream&& out, sim::MCShower const& shower,
    std::string indent = "", bool bIndentFirst = true
    ) const;


    private:

  art::InputTag fInputShowers;    ///< name of MCShower's data product
  std::string fOutputCategory;    ///< name of the stream for output
  unsigned int fDaughtersPerLine; ///< number of daughter IDs printed per line

}; // class sim::DumpMCShowers


//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
sim::DumpMCShowers::DumpMCShowers(Parameters const& config)
  : EDAnalyzer(config)
  , fInputShowers(config().InputShowers())
  , fOutputCategory(config().OutputCategory())
  , fDaughtersPerLine(config().DaughtersPerLine())
{}

//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpMCShowers::DumpMCShower(
  Stream&& out, sim::MCShower const& shower,
  std::string indent /* = "" */, bool bIndentFirst /* = true */
) const {
  if (bIndentFirst) out << indent;
  out
    << "from GEANT track ID=" << shower.TrackID()
    << " PDG ID=" << shower.PdgCode()
    << " from " << OriginDescription(shower.Origin())
    << " via '" << shower.Process() << "'";
  out << "\n" << indent
    << "  starting at ";
  ::PrintMCStep(out, shower.Start());
  out << "\n" << indent
    << "  ending at ";
  ::PrintMCStep(out, shower.End());

  TVector3 const& startDir = shower.StartDir();
  out << "\n" << indent
    << "pointing toward ("
    << startDir.X() << ", " << startDir.Y() << ", " << startDir.Z() << ") cm";
  std::vector<double> const& charges = shower.Charge();
  std::vector<double> const& dQdx = shower.dQdx();
  size_t const nQPlanes = dQdx.size(), nChPlanes = charges.size();
  size_t const nPlanes = std::max(nQPlanes, nChPlanes);
  out << "\n" << indent;
  if (nPlanes > 0) {
    out
      << "dE/dx=" << shower.dEdx() << " MeV/cm and dQ/dx (charge) on "
      << nPlanes << " planes:";
    for (size_t iPlane = 0; iPlane < nPlanes; ++iPlane) {
      out << " [#" << iPlane << "] ";
      if (iPlane < dQdx.size()) out << dQdx[iPlane];
      else                      out << "<N/A>";
      if (iPlane < charges.size()) out << " (" << charges[iPlane] << ")";
      else                      out << "<N/A>";
    } // for plane
  }
  else out << "no energy or charge information available";

  std::vector<unsigned int> const& daughters = shower.DaughterTrackID();
  out << "\n" << indent
    << "combined energy deposition information: ";
  ::PrintMCStep(out, shower.DetProfile());
  out << "\n" << indent
    << daughters.size() << " daughters, ID:";
  for (size_t i = 0; i < daughters.size(); ++i) {
    if ((i % fDaughtersPerLine) == 0) out << "\n" << indent << "    ";
    out << " " << std::setw(8) << daughters[i];
  } // for

  out << "\n" << indent
    << "mother ID=" << shower.MotherTrackID()
    << " PDG ID=" << shower.MotherPdgCode()
    << " via '" << shower.MotherProcess() << "'";
  out << "\n" << indent
    << "  starting at ";
  ::PrintMCStep(out, shower.MotherStart());
  out << "\n" << indent
    << "  ending at ";
  ::PrintMCStep(out, shower.MotherEnd());

  out << "\n" << indent
    << "ancestor ID=" << shower.AncestorTrackID()
    << " PDG ID=" << shower.AncestorPdgCode()
    << " via '" << shower.AncestorProcess() << "'";
  out << "\n" << indent
    << "  starting at ";
  ::PrintMCStep(out, shower.AncestorStart());
  out << "\n" << indent
    << "  ending at ";
  ::PrintMCStep(out, shower.AncestorEnd());

} // sim::DumpMCShowers::DumpMCShower()


//------------------------------------------------------------------------------
void sim::DumpMCShowers::analyze(art::Event const& event) {

  // get the particles from the event
  auto const& Showers
    = *(event.getValidHandle<std::vector<sim::MCShower>>(fInputShowers));

  mf::LogVerbatim(fOutputCategory)
    << "Event " << event.id() << ": data product '"
    << fInputShowers.encode() << "' contains "
    << Showers.size() << " MCShower objects";

  unsigned int iShower = 0;
  mf::LogVerbatim log(fOutputCategory);
  for (sim::MCShower const& shower: Showers) {

    // a bit of a header
    log << "\n[#" << (iShower++) << "] ";
    DumpMCShower(log, shower, "  ", false);

  } // for
  log << "\n";

} // sim::DumpMCShowers::analyze()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpMCShowers)

//------------------------------------------------------------------------------
