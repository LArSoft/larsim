/**
 * @file   DumpSimChannels_module.cc
 * @brief  Module dumping SimChannels information on screen
 * @date   March 30, 2016
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 */


// nutools libraries
#include "lardataobj/Simulation/SimChannel.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


namespace sim {
  class DumpSimChannels;
} // namespace sim

namespace {
  using namespace fhicl;

  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> InputSimChannels {
      Name("InputSimChannels"),
      Comment("data product with the SimChannels to be dumped")
      };

    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the output stream (managed by the message facility)"),
      "DumpSimChannels" /* default value */
      };

  }; // struct Config


} // local namespace


class sim::DumpSimChannels: public art::EDAnalyzer {
    public:
  // type to enable module parameters description by art
  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Configuration-checking constructor
  explicit DumpSimChannels(Parameters const& config);

  // Plugins should not be copied or assigned.
  DumpSimChannels(DumpSimChannels const&) = delete;
  DumpSimChannels(DumpSimChannels &&) = delete;
  DumpSimChannels& operator = (DumpSimChannels const&) = delete;
  DumpSimChannels& operator = (DumpSimChannels &&) = delete;


  // Operates on the event
  void analyze(art::Event const& event) override;


  /**
   * @brief Dumps the content of the specified SimChannel in the output stream
   * @tparam Stream the type of output stream
   * @param out the output stream
   * @param simchannel the SimChannel to be dumped
   * @param indent base indentation string (default: none)
   * @param bIndentFirst if first output line should be indented (default: yes)
   *
   * The indent string is prepended to every line of output, with the possible
   * exception of the first one, in case bIndentFirst is true.
   *
   * The output starts on the current line, and the last line is NOT broken.
   */
  template <typename Stream>
  void DumpSimChannel(
    Stream&& out, sim::SimChannel const& simchannel,
    std::string indent = "", bool bIndentFirst = true
    ) const;


    private:

  art::InputTag fInputChannels; ///< name of SimChannel's data product
  std::string fOutputCategory; ///< name of the stream for output

}; // class sim::DumpSimChannels


//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
sim::DumpSimChannels::DumpSimChannels(Parameters const& config)
  : EDAnalyzer(config)
  , fInputChannels(config().InputSimChannels())
  , fOutputCategory(config().OutputCategory())
{}


//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpSimChannels::DumpSimChannel(
  Stream&& out, sim::SimChannel const& channel,
  std::string indent /* = "" */, bool bIndentFirst /* = true */
) const {
  if (bIndentFirst) out << indent;
  channel.Dump(out, indent);
} // sim::DumpSimChannels::DumpSimChannels()


//------------------------------------------------------------------------------
void sim::DumpSimChannels::analyze(art::Event const& event) {

  // get the particles from the event
  auto const& SimChannels
    = *(event.getValidHandle<std::vector<sim::SimChannel>>(fInputChannels));

  mf::LogVerbatim(fOutputCategory) << "Event " << event.id()
    << " : data product '" << fInputChannels.encode() << "' contains "
    << SimChannels.size() << " SimChannels";

  unsigned int iSimChannel = 0;
  for (sim::SimChannel const& simChannel: SimChannels) {

    // a bit of a header
    mf::LogVerbatim log(fOutputCategory);
    log << "[#" << (iSimChannel++) << "] ";
    DumpSimChannel(log, simChannel, "  ", false);

  } // for
  mf::LogVerbatim(fOutputCategory) << "\n";

} // sim::DumpSimChannels::analyze()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpSimChannels)

//------------------------------------------------------------------------------
