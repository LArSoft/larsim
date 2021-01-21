/**
 * @file   larsim/EventGenerator/POTaccumulator_module.cc
 * @brief  Module summing all POT from the input files.
 * @date   April 29, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

// LArSoft libraries
#include "larcoreobj/SummaryData/POTSummary.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/SummedValue.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "range/v3/view.hpp"

// C/C++ standard libraries
#include <map>
#include <string>

// -----------------------------------------------------------------------------
namespace sim {
  class POTaccumulator;
}

/**
 * @brief  Prints on console the total Protons On Target from the input subruns.
 *
 * This module collects information from each of the subrun in the input files
 * and prints the total of the protons on target. It also optionally prints the
 * total for each observed run.
 * If a subrun is met more than once, the information from all subrun fragments
 * are added together (i.e. it is assumed that summary information is
 * complementary rather than duplicate).
 *
 * The output can be printed to the console or a file via the appropriate
 * message facility configuration.
 *
 * Two output streams are used: the one for the run summary and the one for
 * the total summary. They may coincide.
 *
 *
 * Input
 * ------
 *
 * The module reads information from objects of type `sumdata::POTSummary`
 * stored in each _subrun_.
 *
 *
 * Configuration
 * --------------
 *
 * The output streams can be managed by configuring message facility to deal
 * with the relevant category names. All messages are sent with the INFO level.
 *
 * * `SummaryTag` (input tag; default: `generator`): data product (subrun level)
 *   with the summary information;
 * * `SummaryCategory` (string; default: `POTaccumulator`): the name of the
 *   output category the summary is sent to;
 * * `RunSummaryCategory` (string; default: disabled): if specified, a summary
 *   POT is printed for each run; the summary is sent to the output stream
 *   specified by the value of this parameter.
 *
 *
 */
class sim::POTaccumulator : public art::EDAnalyzer {
public:
  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> SummaryTag{
      Name("SummaryTag"),
      Comment("data product (subrun level) with the summary information"),
      "generator"};

    fhicl::Atom<std::string> SummaryCategory{
      Name("SummaryCategory"),
      Comment("name of the output category the summary is sent to"),
      "POTaccumulator" // default value
    };

    fhicl::Atom<std::string> RunSummaryCategory{
      Name("RunSummaryCategory"),
      Comment("name of the output category the summary is sent to"),
      "" // default value
    };

  }; // struct Config

  /// Type to enable module parameters description by _art_.
  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Configuration-checking constructor.
  explicit POTaccumulator(Parameters const& config);

  // Plugins should not be copied or assigned.
  POTaccumulator(POTaccumulator const&) = delete;
  POTaccumulator(POTaccumulator&&) = delete;
  POTaccumulator& operator=(POTaccumulator const&) = delete;
  POTaccumulator& operator=(POTaccumulator&&) = delete;

  // Nothing to be done at event level.
  virtual void
  analyze(art::Event const& event) override
  {}

  /// Collects information from each subrun.
  virtual void endSubRun(art::SubRun const& subRun) override;

  /// Prints the general summary.
  virtual void endJob() override;

private:
  // -- BEGIN -- Configuration variables ---------------------------------------

  art::InputTag fPOTtag;              ///< Name of `sumdata::POTSummary` data product.
  std::string fSummaryOutputCategory; ///< Name of the main stream for output.
  std::string fRunOutputCategory;     ///< Name of the run stream for output.

  // -- END -- Configuration variables -----------------------------------------

  // -- BEGIN -- Internal cache variables --------------------------------------

  /// Count of subrun fragments with POT information.
  std::map<art::SubRunID, unsigned int> fPresentSubrunFragments;

  /// Count of subrun fragments without POT information.
  std::map<art::SubRunID, unsigned int> fMissingSubrunFragments;

  /// Partial count of POT in the run, per run.
  std::map<art::RunID, art::SummedValue<sumdata::POTSummary>> fRunPOT;

  // -- END -- Internal cache variables ----------------------------------------

  /// Prints the list of subruns with partial or missing POT information.
  void printMissingSubrunList() const;

  /// Prints the list of POT per run.
  void printRunSummary() const;

  /// Prints the total POT summary `totalPOT`.
  void printSummary(sumdata::POTSummary const& totalPOT) const;

  /// Converts the information from `POT` in a compact string.
  static std::string to_string(sumdata::POTSummary const& POT);

}; // class sim::POTaccumulator

//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
sim::POTaccumulator::POTaccumulator(Parameters const& config)
  : EDAnalyzer(config)
  , fPOTtag(config().SummaryTag())
  , fSummaryOutputCategory(config().SummaryCategory())
  , fRunOutputCategory(config().RunSummaryCategory())
{}

//------------------------------------------------------------------------------
void
sim::POTaccumulator::endSubRun(art::SubRun const& subRun)
{

  auto const& ID = subRun.id();

  //
  // get the information from the subrun and update the subrun counts
  //
  art::Handle<sumdata::POTSummary> summaryHandle;
  if (!subRun.getByLabel(fPOTtag, summaryHandle)) {
    ++fMissingSubrunFragments[ID];
    mf::LogDebug(fSummaryOutputCategory)
      << "Fragment of subrun " << ID << " has no '" << fPOTtag.encode() << "' POT summary.";
    return;
  }

  ++fPresentSubrunFragments[ID];

  //
  // accumulate the information by run
  //
  sumdata::POTSummary const& subRunPOT = *summaryHandle;

  fRunPOT[ID.runID()].update(summaryHandle);
  MF_LOG_TRACE(fSummaryOutputCategory)
    << "Fragment #" << fPresentSubrunFragments[ID] << " of subrun " << ID << ": "
    << sim::POTaccumulator::to_string(subRunPOT);

} // sim::POTaccumulator::endSubRun()

//------------------------------------------------------------------------------
void
sim::POTaccumulator::endJob()
{

  //
  // print the run summary
  //

  if (!fRunOutputCategory.empty()) {

    if (!fMissingSubrunFragments.empty()) printMissingSubrunList();

    printRunSummary();

  } // if

  //
  // print the total summary
  //

  // here we skip _art_ aggregation mechanism
  // because it can't handle multiple runs
  sumdata::POTSummary totalPOT;
  for (auto const& POT : fRunPOT | ranges::views::values)
    totalPOT.aggregate(POT.value());

  printSummary(totalPOT);

} // sim::POTaccumulator::endJob()

//------------------------------------------------------------------------------
void
sim::POTaccumulator::printMissingSubrunList() const
{

  //
  // missing fragments information
  //
  mf::LogVerbatim log{fRunOutputCategory};
  log << size(fMissingSubrunFragments) << " subruns lack POT information:";

  auto const fend = fPresentSubrunFragments.cend();

  for (auto const& [id, nMissing] : fMissingSubrunFragments) {

    // add to the count of fragments the ones which we have actually found
    unsigned int nFragments = nMissing;
    auto const iFound = fPresentSubrunFragments.find(id);
    if (iFound != fend) nFragments += iFound->second;

    log << "\n" << id << ": " << nMissing << " / " << nFragments << " \"fragments\"";

  } // for

} // sim::POTaccumulator::printMissingSubrunList()

//------------------------------------------------------------------------------
void
sim::POTaccumulator::printRunSummary() const
{

  // count subruns in run
  std::map<art::RunID, unsigned int> subrunCount;
  for (art::SubRunID const& ID : fPresentSubrunFragments | ranges::views::keys)
    ++subrunCount[ID.runID()];

  mf::LogVerbatim log{fRunOutputCategory};
  log << "POT from " << size(fRunPOT) << " runs:";
  for (auto const& [id, POT] : fRunPOT) {
    log << "\n " << id << " (" << subrunCount[id]
        << " subruns): " << sim::POTaccumulator::to_string(POT.value());
  } // for

} // sim::POTaccumulator::printRunSummary()

//------------------------------------------------------------------------------
void
sim::POTaccumulator::printSummary(sumdata::POTSummary const& totalPOT) const
{

  // aggregate all run summaries
  mf::LogVerbatim{fSummaryOutputCategory}
    << "Aggregated POT from " << fRunPOT.size() << " runs (" << fPresentSubrunFragments.size()
    << " subruns): " << sim::POTaccumulator::to_string(totalPOT);

} // sim::POTaccumulator::printSummary()

//------------------------------------------------------------------------------
std::string
sim::POTaccumulator::to_string(sumdata::POTSummary const& POT)
{
  using namespace std::string_literals;
  return std::to_string(POT.totgoodpot) + " good POT ( "s + std::to_string(POT.goodspills) +
         " spills); total: " + std::to_string(POT.totpot) + " POT ( "s +
         std::to_string(POT.totspills) + " spills)";
} // sim::POTaccumulator::to_string(sumdata::POTSummary)

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::POTaccumulator)

//------------------------------------------------------------------------------
