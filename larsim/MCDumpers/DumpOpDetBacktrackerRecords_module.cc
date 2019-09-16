/**
 * @file   DumpOpDetBacktrackerRecords_module.cc
 * @brief  Module dumping OpDetBacktrackerRecords information on screen
 * @date   August 8, 2016
 * @author Lucas Mendes Santos (lmendes@ifi.unicamp.br)
 *
 */


// nutools libraries
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


namespace sim {
  class DumpOpDetBacktrackerRecords;
} // namespace sim

namespace {
  using namespace fhicl;

  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> InputOpDetBacktrackerRecord {
      Name("InputOpDetBacktrackerRecord"),
      Comment("data product with the OpDetBacktrackerRecord to be dumped")
      };

    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the output stream (managed by the message facility)"),
      "DumpOpDetBacktrackerRecords" /* default value */
      };

  }; // struct Config


} // local namespace


class sim::DumpOpDetBacktrackerRecords: public art::EDAnalyzer {
    public:
  // type to enable module parameters description by art
  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Configuration-checking constructor
  explicit DumpOpDetBacktrackerRecords(Parameters const& config);

  // Plugins should not be copied or assigned.
  DumpOpDetBacktrackerRecords(DumpOpDetBacktrackerRecords const&) = delete;
  DumpOpDetBacktrackerRecords(DumpOpDetBacktrackerRecords &&) = delete;
  DumpOpDetBacktrackerRecords& operator = (DumpOpDetBacktrackerRecords const&) = delete;
  DumpOpDetBacktrackerRecords& operator = (DumpOpDetBacktrackerRecords &&) = delete;


  // Operates on the event
  void analyze(art::Event const& event) override;


  /**
   * @brief Dumps the content of the specified OpDetBacktrackerRecord in the output stream
   * @tparam Stream the type of output stream
   * @param out the output stream
   * @param simchannel the OpDetBacktrackerRecord to be dumped
   * @param indent base indentation string (default: none)
   * @param bIndentFirst if first output line should be indented (default: yes)
   *
   * The indent string is prepended to every line of output, with the possible
   * exception of the first one, in case bIndentFirst is true.
   *
   * The output starts on the current line, and the last line is NOT broken.
   */
  template <typename Stream>
  void DumpOpDetBacktrackerRecord(
    Stream&& out, sim::OpDetBacktrackerRecord const& simchannel,
    std::string indent = "", bool bIndentFirst = true
    ) const;


    private:

  art::InputTag fInputChannels; ///< name of OpDetBacktrackerRecord's data product
  std::string fOutputCategory; ///< name of the stream for output

}; // class sim::DumpOpDetBacktrackerRecords


//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
sim::DumpOpDetBacktrackerRecords::DumpOpDetBacktrackerRecords(Parameters const& config)
  : EDAnalyzer(config)
  , fInputChannels(config().InputOpDetBacktrackerRecord())
  , fOutputCategory(config().OutputCategory())
{}


//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpOpDetBacktrackerRecords::DumpOpDetBacktrackerRecord(
  Stream&& out, sim::OpDetBacktrackerRecord const& channel,
  std::string indent /* = "" */, bool bIndentFirst /* = true */
) const {
  if (bIndentFirst) out << indent;
  channel.Dump(out, indent);
} // sim::DumpOpDetBacktrackerRecords::DumpOpDetBacktrackerRecords()


//------------------------------------------------------------------------------
void sim::DumpOpDetBacktrackerRecords::analyze(art::Event const& event) {

  // get the particles from the event
  auto const& OpDetBacktrackerRecord
    = *(event.getValidHandle<std::vector<sim::OpDetBacktrackerRecord>>(fInputChannels));

  mf::LogVerbatim(fOutputCategory) << "Event " << event.id()
    << " : data product '" << fInputChannels.encode() << "' contains "
    << OpDetBacktrackerRecord.size() << " OpDetBacktrackerRecord";

  unsigned int iOpDetBacktrackerRecord = 0;
  for (sim::OpDetBacktrackerRecord const& simChannel: OpDetBacktrackerRecord) {

    // a bit of a header
    mf::LogVerbatim log(fOutputCategory);
    log << "[#" << (iOpDetBacktrackerRecord++) << "] ";
    DumpOpDetBacktrackerRecord(log, simChannel, "  ", false);

  } // for
  mf::LogVerbatim(fOutputCategory) << "\n";

} // sim::DumpOpDetBacktrackerRecords::analyze()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpOpDetBacktrackerRecords)

//------------------------------------------------------------------------------
