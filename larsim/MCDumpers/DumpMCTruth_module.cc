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
#include "MCDumperUtils.h"

// nutools libraries
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


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
  
  
  /**
   * @brief Dumps the content of the specified MC truth in the output stream.
   * @tparam Stream the type of output stream
   * @param out the output stream
   * @param truth the truth information to be dumped
   * @param indent base indentation string (default: none)
   * @param firstIndent string to be used for indentation of the first line
   * 
   * The indent string is prepended to every line of output, with the possible
   * exception of the first one, in case bIndentFirst is true.
   * 
   * The output starts on the current line, and the last line is NOT broken.
   */
  template <typename Stream>
  void DumpTruth(
    Stream&& out, simb::MCTruth const& truth,
    std::string indent = "", std::string firstIndent = ""
    ) const;
  
  /// Returns the name of the product in the form `"module_instance_process"`.
  template <typename Handle>
  static std::string productName(Handle const& handle);
  
  
    private:
  
  std::vector<art::InputTag> fInputTruth; ///< Name of MCTruth data products.
  std::string fOutputCategory; ///< Name of the stream for output.
  bool bAllTruth = false; ///< Whether to process all MCTruth collections.
  
  unsigned int fPointsPerLine; ///< trajectory points per output line
  
  
  /**
   * @brief Dumps the content of the specified particle in the output stream.
   * @tparam Stream the type of output stream
   * @param out the output stream
   * @param particle the particle to be dumped
   * @param indent base indentation string (default: none)
   * @param firstIndent string to be used for indentation of the first line
   * 
   * The indent string is prepended to every line of output, with the possible
   * exception of the first one, in case bIndentFirst is true.
   * 
   * The output starts on the current line, and the last line is NOT broken.
   */
  template <typename Stream>
  void DumpMCParticle(
    Stream&& out, simb::MCParticle const& particle,
    std::string indent = "", std::string firstIndent = ""
    ) const;
  
  /**
   * @brief Dumps the content of the specified neutrino in the output stream.
   * @tparam Stream the type of output stream
   * @param out the output stream
   * @param neutrino the neutrino to be dumped
   * @param indent base indentation string (default: none)
   * @param firstIndent string to be used for indentation of the first line
   * 
   * The indent string is prepended to every line of output, with the possible
   * exception of the first one, in case bIndentFirst is true.
   * 
   * The output starts on the current line, and the last line is NOT broken.
   */
  template <typename Stream>
  void DumpMCNeutrino(
    Stream&& out, simb::MCNeutrino const& neutrino,
    std::string indent = "", std::string firstIndent = ""
    ) const;
  
  
}; // class sim::DumpMCTruth


//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
namespace {
  
  template <typename Stream>
  Stream& operator<< (Stream&& out, TLorentzVector const& v) {
    out
      << "(" << v.X() << ", " << v.Y() << ", " << v.Z() << "; " << v.T() << ")";
    return out;
  } // operator<< (TLorentzVector)
  
} // local namespace


//------------------------------------------------------------------------------
sim::DumpMCTruth::DumpMCTruth(Parameters const& config)
  : EDAnalyzer(config)
  , fInputTruth()
  , fOutputCategory(config().OutputCategory())
  , bAllTruth(!config().InputTruth(fInputTruth)) // true if InputTruth omitted
  , fPointsPerLine(config().PointsPerLine())
  {}

//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpMCTruth::DumpMCParticle(
  Stream&& out, simb::MCParticle const& particle,
  std::string indent /* = "" */, std::string firstIndent /* = "" */
) const {
  out << firstIndent
    << ParticleName(particle.PdgCode())
    << "  (ID=" << particle.TrackId() << ")"
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
    out << "\n" << indent << "created at " << start << " cm by "
      << (particle.Process().empty()? "magics": particle.Process())
      << " from ID=" << particle.Mother()
      << " with momentum " << start_mom << " GeV/c";
    out << "\n" << indent;
    if (nPoints == 1) {
      out << "still alive!";
    }
    else {
      TLorentzVector const& stop = particle.EndPosition();
      TLorentzVector const& stop_mom = particle.EndMomentum();
      out << "stops at " << stop << " cm with momentum " << stop_mom
        << " GeV/c by "
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
    
    out << "\n" << indent << "comes with a trajectory";
    if (nPoints == 1) {
      out << " of a single point";
      if (fPointsPerLine > 0) out << ": " << traj[0].first << " cm";
    }
    else {
      out
        << " " << traj.TotalLength() << " cm long in " << nPoints << " points";
      if (fPointsPerLine > 0) {
        out << ":";
        unsigned int page = 0;
        for (auto const& pair: traj) {
          if (page-- == 0) {
            out << "\n" << indent << "  ";
            page = fPointsPerLine - 1;
          }
          else out << " -- ";
          
          out << pair.first; // position
        } // for trajectory points
      } // if printing points
    }
  } // if has points
  
} // sim::DumpMCTruth::DumpMCParticle()


//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpMCTruth::DumpMCNeutrino(
  Stream&& out, simb::MCNeutrino const& nu,
  std::string indent /* = "" */, std::string firstIndent /* = "" */
) const {
  out << firstIndent
    << "from " << TruthCCNCname(nu.CCNC())
    << ", " << TruthInteractionTypeName(nu.InteractionType())
    << ", mode: " << nu.Mode()
    << ", target: " << nu.Target()
    << ", hit nucleon: " << nu.HitNuc()
    << ", hit quark: " << nu.HitQuark()
    << '\n' << indent
    << "x=" << nu.X() << " y=" << nu.Y() << " w=" << nu.W()
    << " Q^2=" << nu.QSqr() << " GeV^2; theta=" << nu.Theta()
    << " rad pT=" << nu.Pt() << " GeV/c";
  out << '\n' << indent << "neutrino: ";
  DumpMCParticle(std::forward<Stream>(out), nu.Nu(), indent + "  ");
  out << '\n' << indent << "outgoing lepton: ";
  DumpMCParticle(std::forward<Stream>(out), nu.Lepton(), indent + "  ");
  
} // sim::DumpMCTruth::DumpMCNeutrino()


//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpMCTruth::DumpTruth(
  Stream&& out, simb::MCTruth const& truth,
  std::string indent /* = "" */, std::string firstIndent /* = "" */
) const {
  unsigned int const nParticles = truth.NParticles();
  out << firstIndent
    << nParticles << " particles from "
    << TruthOriginName(truth.Origin());
  if (truth.NeutrinoSet()) {
    out << '\n' << indent << "neutrino information: ";
    DumpMCNeutrino
      (std::forward<Stream>(out), truth.GetNeutrino(), indent + "    ");
  }
  for (unsigned int i = 0; i < nParticles; ++i) {
    out << '\n' << indent << "[#" << i << "] ";
    DumpMCParticle(std::forward<Stream>
      (out), truth.GetParticle(i), indent + "  ", "");
  } // for all particles
  
} // sim::DumpMCTruth::DumpTruth()


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
      DumpTruth(log, truth, "  ", "");
      
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
