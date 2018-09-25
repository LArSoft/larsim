/**
 * @file   larsim/MCDumpers/DumpGTruth_module.cc
 * @brief  Module dumping GTruth information from GENIE on screen.
 * @date   November 2, 2017
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * 
 * This modules complements `DumpMCTruth` module for GENIE output.
 * 
 */

// LArSoft libraries
#include "lardataalg/MCDumpers/MCDumpers.h" // sim::dump namespace

// nutools libraries
#include "nusimdata/SimulationBase/GTruth.h"

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

// C/C++ standard libraries
#include <vector>
#include <string>
#include <iterator> // std::back_inserter()
#include <algorithm> // std::copy()
#include <utility> // std::forward()

namespace sim {
  class DumpGTruth;
}

class sim::DumpGTruth: public art::EDAnalyzer {
    public:
  
  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::OptionalSequence<art::InputTag> InputTruth {
      Name("InputTruth"),
      Comment("data product with the collection of GENIE truth to be dumped")
      };
    
    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the output stream (managed by the message facility)"),
      "DumpGTruth" /* default value */
      };
    
    fhicl::Atom<bool> AllowNoTruth {
      Name("AllowNoTruth"),
      Comment("when InputTruth is empty, allow for no truth to be found"),
      false /* default value */
      };
    
  }; // struct Config
  
  
  /// Type to enable module parameters description by _art_.
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  /// Configuration-checking constructor.
  explicit DumpGTruth(Parameters const& config);
  
  // Plugins should not be copied or assigned.
  DumpGTruth(DumpGTruth const&) = delete;
  DumpGTruth(DumpGTruth &&) = delete;
  DumpGTruth& operator = (DumpGTruth const&) = delete;
  DumpGTruth& operator = (DumpGTruth &&) = delete;
  
  
  // Operates on the event
  virtual void analyze(art::Event const& event) override;
  
  
  /// Returns the name of the product in the form `"module_instance_process"`.
  template <typename Handle>
  static std::string productName(Handle const& handle);
  
    private:
  
  std::vector<art::InputTag> fInputTruth; ///< Name of `GTruth` data products.
  std::string fOutputCategory; ///< Name of the stream for output.
  bool bAllTruth = false; ///< Whether to process all `GTruth` collections.
  bool bAllowNoTruth = false; ///< Whether to forgive when no truth is present.
  
}; // class sim::DumpGTruth


//------------------------------------------------------------------------------
//---  module implementation
//------------------------------------------------------------------------------
sim::DumpGTruth::DumpGTruth(Parameters const& config)
  : EDAnalyzer(config)
  , fInputTruth()
  , fOutputCategory(config().OutputCategory())
  , bAllTruth(!config().InputTruth(fInputTruth)) // true if InputTruth omitted
  , bAllowNoTruth(config().AllowNoTruth())
{
  if (!bAllTruth && bAllowNoTruth) {
    throw art::Exception(art::errors::Configuration)
      << "'AllowNoTruth' is only allowed if no 'InputTruth' is specified.\n";
  }
}


//------------------------------------------------------------------------------
void sim::DumpGTruth::analyze(art::Event const& event) {
  
  //
  // prepare the data products to be dumped
  //
  struct ProductInfo_t {
    using Thruths_t = std::vector<simb::GTruth>;
    Thruths_t const* truths;
    std::string name;
    
    ProductInfo_t(art::Handle<Thruths_t> const& handle)
      : truths(handle.provenance()->isPresent()? handle.product(): nullptr)
      , name(sim::DumpGTruth::productName(handle))
      {}
    ProductInfo_t(art::ValidHandle<Thruths_t> const& handle)
      : truths(handle.product()), name(sim::DumpGTruth::productName(handle))
      {}
    
  }; // ProductInfo_t
  
  std::vector<ProductInfo_t> AllTruths;
  if (bAllTruth) {
    std::vector<art::Handle<std::vector<simb::GTruth>>> handles;
    event.getManyByType(handles);
    std::copy(handles.begin(), handles.end(), std::back_inserter(AllTruths));
  }
  else {
    for (auto const& inputTag: fInputTruth) {
      AllTruths.emplace_back
        (event.getValidHandle<std::vector<simb::GTruth>>(inputTag));
    } // for
  }
  
  //
  // sanity check
  //
  if (AllTruths.empty() && !bAllowNoTruth) {
    throw art::Exception(art::errors::ProductNotFound)
      << "No GENIE truth found to be dumped!\n";
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
      << " contains " << nTruths << " GENIE truth blocks in "
      << AllTruths.size() << " collections";
  }
  else if (AllTruths.size() == 1) {
    mf::LogVerbatim(fOutputCategory) << "Event " << event.id();
  }
  else {
    mf::LogVerbatim(fOutputCategory) << "Dumping " << nTruths
      << " GENIE truth blocks from " << AllTruths.size()
      << " collections in event " << event.id();
  }
  
  //
  // dump data product by data product
  //
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
    // dump each GENIE truth in the data product
    //
    unsigned int iTruth = 0;
    for (auto const& truth: *truths) {
      
      mf::LogVerbatim log (fOutputCategory);
      
      if (truths->size() > 1) log << "(#" << iTruth << ") ";
      sim::dump::DumpGTruth(log, truth, "  ", "");
      
      //
      // update counters
      //
      ++iTruth;
      
    } // for each truth in data product
    
  } // for truth data products
  
  //
  // all done
  //
  
} // sim::DumpGTruth::analyze()


//------------------------------------------------------------------------------
template <typename Handle>
std::string sim::DumpGTruth::productName(Handle const& handle) {
  auto const* prov = handle.provenance();
  return prov->moduleLabel()
    + '_' + prov->productInstanceName()
    + '_' + prov->processName()
    ;
} // sim::DumpGTruth::productName()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpGTruth)

//------------------------------------------------------------------------------
