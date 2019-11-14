/**
 * @file   DumpMCTruth_module.cc
 * @brief  Module dumping MCTruth information on screen.
 * @date   February 22nd, 2017
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 * Dedication: to the memory of DumpMCTruth Sr., whose untimely departure by a
 * rogue `rm -R` deprived us of a useful tool, and the author of a good deal of
 * good karma.
 *
 */

// LArSoft libraries
#include "lardataalg/MCDumpers/MCDumpers.h" // sim::dump namespace

// nusimdata libraries
#include "nusimdata/SimulationBase/MCTruth.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace sim {
  class DumpMCTruth;
}


class sim::DumpMCTruth: public art::EDAnalyzer {
    public:

  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::OptionalSequence<art::InputTag> InputTruth {
      Name("InputTruth"),
      Comment("data product with the collection of MC truth to be dumped")
      };

    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the output stream (managed by the message facility)"),
      "DumpMCTruth" /* default value */
      };

    fhicl::Atom<unsigned int> PointsPerLine {
      Name("PointsPerLine"),
      Comment("trajectory points printed per line (default: 2; 0 = skip them)"),
      2 /* default value */
      };

  }; // struct Config

  /// Type to enable module parameters description by _art_.
  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Configuration-checking constructor.
  explicit DumpMCTruth(Parameters const& config);

  // Plugins should not be copied or assigned.
  DumpMCTruth(DumpMCTruth const&) = delete;
  DumpMCTruth(DumpMCTruth &&) = delete;
  DumpMCTruth& operator = (DumpMCTruth const&) = delete;
  DumpMCTruth& operator = (DumpMCTruth &&) = delete;


  // Operates on the event
  void analyze(art::Event const& event) override;


  /// Returns the name of the product in the form `"module_instance_process"`.
  template <typename Handle>
  static std::string productName(Handle const& handle);


    private:

  std::vector<art::InputTag> fInputTruth; ///< Name of MCTruth data products.
  std::string fOutputCategory; ///< Name of the stream for output.
  bool bAllTruth = false; ///< Whether to process all MCTruth collections.

  unsigned int fPointsPerLine; ///< trajectory points per output line


}; // class sim::DumpMCTruth


//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
sim::DumpMCTruth::DumpMCTruth(Parameters const& config)
  : EDAnalyzer(config)
  , fInputTruth()
  , fOutputCategory(config().OutputCategory())
  , bAllTruth(!config().InputTruth(fInputTruth)) // true if InputTruth omitted
  , fPointsPerLine(config().PointsPerLine())
  {}


//------------------------------------------------------------------------------
void sim::DumpMCTruth::analyze(art::Event const& event) {

  //
  // prepare the data products to be dumped
  //
  struct ProductInfo_t {
    using Thruths_t = std::vector<simb::MCTruth>;
    Thruths_t const* truths;
    std::string name;

    ProductInfo_t(art::Handle<Thruths_t> const& handle)
      : truths(handle.provenance()->isPresent()? handle.product(): nullptr)
      , name(sim::DumpMCTruth::productName(handle))
      {}
    ProductInfo_t(art::ValidHandle<Thruths_t> const& handle)
      : truths(handle.product()), name(sim::DumpMCTruth::productName(handle))
      {}

  }; // ProductInfo_t

  std::vector<ProductInfo_t> AllTruths;
  if (bAllTruth) {
    std::vector<art::Handle<std::vector<simb::MCTruth>>> handles;
    event.getManyByType(handles);
    std::copy(handles.begin(), handles.end(), std::back_inserter(AllTruths));
  }
  else {
    for (auto const& inputTag: fInputTruth) {
      AllTruths.emplace_back
        (event.getValidHandle<std::vector<simb::MCTruth>>(inputTag));
    } // for
  }

  //
  // sanity check
  //
  if (AllTruths.empty()) {
    throw art::Exception(art::errors::ProductNotFound)
      << "No MC truth found to be dumped!\n";
  }

  //
  // print an introduction
  //
  unsigned int const nTruths = std::accumulate(
    AllTruths.begin(), AllTruths.end(), 0U,
    [](unsigned int total, auto const& info)
      { return total + (info.truths? info.truths->size(): 0); }
    );

  if (bAllTruth) {
    mf::LogVerbatim(fOutputCategory) << "Event " << event.id()
      << " contains " << nTruths << " MC truth blocks in "
      << AllTruths.size() << " collections";
  }
  else if (AllTruths.size() == 1) {
    mf::LogVerbatim(fOutputCategory) << "Event " << event.id();
  }
  else {
    mf::LogVerbatim(fOutputCategory) << "Dumping " << nTruths
      << " MC truth blocks from " << AllTruths.size()
      << " collections in event " << event.id();
  }

  //
  // dump data product by data product
  //
  unsigned int nParticles = 0, nNeutrinos = 0;
  for (ProductInfo_t const& truths_info: AllTruths) {

    auto const* truths = truths_info.truths;
    std::string productName = truths_info.name;

    if (!truths) {
      mf::LogVerbatim(fOutputCategory)
        << "Data product '" << productName
        << "' has been dropped. No information available.";
    }

    if (AllTruths.size() > 1) {
      mf::LogVerbatim(fOutputCategory)
        << "Data product '" << productName
        << "' contains " << truths->size() << " truth blocks:";
    }
    else if (truths->size() > 1) {
      mf::LogVerbatim(fOutputCategory)
        << truths->size() << " truth blocks:";
    }

    //
    // dump each MC truth in the data product
    //
    unsigned int iTruth = 0;
    for (auto const& truth: *truths) {

      mf::LogVerbatim log (fOutputCategory);

      if (truths->size() > 1) log << "(#" << iTruth << ") ";
      sim::dump::DumpMCTruth(log, truth, "  ", "");

      //
      // update counters
      //
      ++iTruth;
      nParticles += truth.NParticles();
      if (truth.NeutrinoSet()) ++nNeutrinos;

    } // for each truth in data product

  } // for truth data products

  //
  // all done
  //
  mf::LogVerbatim(fOutputCategory) << nNeutrinos
    << " neutrinos generated, " << nParticles
    << " generated particles to be simulated downstream.";

} // sim::DumpMCTruth::analyze()


//------------------------------------------------------------------------------
template <typename Handle>
std::string sim::DumpMCTruth::productName(Handle const& handle) {
  auto const* prov = handle.provenance();
  return prov->moduleLabel()
    + '_' + prov->productInstanceName()
    + '_' + prov->processName()
    ;
} // sim::DumpMCTruth::productName()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpMCTruth)

//------------------------------------------------------------------------------
