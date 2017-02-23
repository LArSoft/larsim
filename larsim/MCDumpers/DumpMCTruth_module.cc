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

// ROOT
#include "TDatabasePDG.h"
#include "TParticlePDG.h"


namespace sim {
  class DumpMCTruth;
} // namespace sim

namespace {
  using namespace fhicl;
  
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
  
  
} // local namespace


class sim::DumpMCTruth: public art::EDAnalyzer {
    public:
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
  
  
  /// Returns a string representing the specified origin.
  static std::string TruthOriginName(simb::Origin_t origin);
  
  /// Returns a string representing the specified process.
  static std::string TruthCCNCname(int ccnc);
  
  /// Returns a string representing the specified interaction type.
  static std::string TruthInteractionTypeName(int type);
  
  /// Returns a string with the name of particle the specified with PDG ID.
  static std::string ParticleName(int pigid);
  
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
    out << "\n" << indent << "created at " << start << " cm by ";
    if (particle.Mother() == 0) out << "the gods";
    else {
      out << (particle.Process().empty()? "magics": particle.Process())
        << " from ID=" << particle.Mother();
    }
    out << " with momentum " << start_mom << " GeV/c";
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
    std::vector<simb::MCTruth> const* truths;
    art::Provenance const* info;
    ProductInfo_t(
      std::vector<simb::MCTruth> const* truths,
      art::Provenance const* info
      )
      : truths(truths), info(info)
      {}
  }; // ProductInfo_t
  std::vector<ProductInfo_t> AllTruths;
  if (bAllTruth) {
    std::vector<art::Handle<std::vector<simb::MCTruth>>> handles;
    event.getManyByType(handles);
    for (auto handle: handles) {
      art::Provenance const* info = handle.provenance();
      AllTruths.emplace_back(
        (info->isPresent()? handle.product(): nullptr),
        info
        );
    } // for handles
  }
  else {
    for (auto const& inputTag: fInputTruth) {
      auto handle = event.getValidHandle<std::vector<simb::MCTruth>>(inputTag);
      AllTruths.emplace_back(
        handle.product(),
        handle.provenance()
        );
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
    auto const* info = truths_info.info;
    std::string productName = info->moduleLabel()
      + '_' + info->productInstanceName()
      + '_' + info->processName()
      ;
    
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
std::string sim::DumpMCTruth::TruthOriginName(simb::Origin_t origin) {
  switch (origin) {
    case simb::kUnknown          : return "unknown origin";
    case simb::kBeamNeutrino     : return "neutrinos from beam";
    case simb::kCosmicRay        : return "cosmic rays";
    case simb::kSuperNovaNeutrino: return "supernova neutrinos";
    case simb::kSingleParticle   : return "single particle";
    default:        return "unsupported (" + std::to_string((int)origin) + ")";
  } // switch
} // sim::DumpMCTruth::TruthOriginName()

//------------------------------------------------------------------------------
std::string sim::DumpMCTruth::TruthCCNCname(int ccnc) {
  switch (ccnc) {
    case simb::kCC: return "charged weak current";
    case simb::kNC: return "neutral weak current";
    default:        return "unsupported (" + std::to_string(ccnc) + ")";
  } // switch
} // sim::DumpMCTruth::TruthCCNCname()


//------------------------------------------------------------------------------
std::string sim::DumpMCTruth::TruthInteractionTypeName(int type) {
  switch (type) {
    case simb::kUnknownInteraction            : return "unknown interaction";
    case simb::kQE                            : return "quasi-elastic scattering";
    case simb::kRes                           : return "resonant scattering";
    case simb::kDIS                           : return "deep inelastic scattering";
    case simb::kCoh                           : return "coherent scattering";
    case simb::kCohElastic                    : return "coherent elastic scattering";
    case simb::kElectronScattering            : return "electron scattering";
    case simb::kIMDAnnihilation               : return "inverse muon decay annihilation";
    case simb::kInverseBetaDecay              : return "inverse beta decay";
    case simb::kGlashowResonance              : return "Glashow resonance";
    case simb::kAMNuGamma                     : return "anomalous neutrino-photon interaction";
    case simb::kMEC                           : return "meson exchange current";
    case simb::kDiffractive                   : return "diffractive";
    case simb::kEM                            : return "electromagnetic";
    case simb::kWeakMix                       : return "weak mixing";
    case simb::kNuanceOffset                  : return "<nuance offset>";
    case simb::kCCQE                          : return "charged current quasi-elastic scattering";
    case simb::kNCQE                          : return "neutral current quasi-elastic scattering";
    case simb::kResCCNuProtonPiPlus           : return "resonant charged current neutrino proton pi+";
    case simb::kResCCNuNeutronPi0             : return "resonant charged current neutrino neutron pi0";
    case simb::kResCCNuNeutronPiPlus          : return "resonant charged current neutrino neutron pi+";
    case simb::kResNCNuProtonPi0              : return "resonant neutral current neutrino proton pi0";
    case simb::kResNCNuProtonPiPlus           : return "resonant neutral current neutrino proton pi+";
    case simb::kResNCNuNeutronPi0             : return "resonant neutral current neutrino neutron pi0";
    case simb::kResNCNuNeutronPiMinus         : return "resonant neutral current neutrino neutron pi-";
    case simb::kResCCNuBarNeutronPiMinus      : return "resonant charged current antineutrino neutron pi-";
    case simb::kResCCNuBarProtonPi0           : return "resonant charged current antineutrino proton pi0";
    case simb::kResCCNuBarProtonPiMinus       : return "resonant charged current antineutrino proton pi-";
    case simb::kResNCNuBarProtonPi0           : return "resonant neutral current antineutrino proton pi0";
    case simb::kResNCNuBarProtonPiPlus        : return "resonant neutral current antineutrino proton pi+";
    case simb::kResNCNuBarNeutronPi0          : return "resonant neutral current antineutrino neutron pi0";
    case simb::kResNCNuBarNeutronPiMinus      : return "resonant neutral current antineutrino neutron pi-";
    case simb::kResCCNuDeltaPlusPiPlus        : return "resonant charged current neutrino Delta+ pi+";
    case simb::kResCCNuDelta2PlusPiMinus      : return "resonant charged current neutrino Delta++ pi-";
    case simb::kResCCNuBarDelta0PiMinus       : return "resonant charged current antineutrino Delta0 pi-";
    case simb::kResCCNuBarDeltaMinusPiPlus    : return "resonant charged current antineutrino Delta- pi+";
    case simb::kResCCNuProtonRhoPlus          : return "resonant charged current neutrino proton rho+";
    case simb::kResCCNuNeutronRhoPlus         : return "resonant charged current neutrino neutron rho+";
    case simb::kResCCNuBarNeutronRhoMinus     : return "resonant charged current antineutrino neutron rho-";
    case simb::kResCCNuBarNeutronRho0         : return "resonant charged current antineutrino neutron rho0";
    case simb::kResCCNuSigmaPlusKaonPlus      : return "resonant charged current neutrino Sigma+ kaon+";
    case simb::kResCCNuSigmaPlusKaon0         : return "resonant charged current neutrino Sigma+ kaon0";
    case simb::kResCCNuBarSigmaMinusKaon0     : return "resonant charged current antineutrino Sigma- kaon0";
    case simb::kResCCNuBarSigma0Kaon0         : return "resonant charged current antineutrino Sigma0 kaon0";
    case simb::kResCCNuProtonEta              : return "resonant charged current neutrino proton eta";
    case simb::kResCCNuBarNeutronEta          : return "resonant charged current antineutrino neutron eta";
    case simb::kResCCNuKaonPlusLambda0        : return "resonant charged current neutrino Kaon+ Lambda0";
    case simb::kResCCNuBarKaon0Lambda0        : return "resonant charged current antineutrino kaon0 Lambda0";
    case simb::kResCCNuProtonPiPlusPiMinus    : return "resonant charged current neutrino proton pi+ pi-";
    case simb::kResCCNuProtonPi0Pi0           : return "resonant charged current neutrino proton pi0 pi0";
    case simb::kResCCNuBarNeutronPiPlusPiMinus: return "resonant charged current antineutrino neutron pi+ pi-";
    case simb::kResCCNuBarNeutronPi0Pi0       : return "resonant charged current antineutrino neutron pi0 pi0";
    case simb::kResCCNuBarProtonPi0Pi0        : return "resonant charged current antineutrino proton pi0 pi0";
    case simb::kCCDIS                         : return "charged current deep inelastic scattering";
    case simb::kNCDIS                         : return "neutral current deep inelastic scattering";
    case simb::kUnUsed1                       : return "unused (1)";
    case simb::kUnUsed2                       : return "unused (2)";
    case simb::kCCQEHyperon                   : return "charged current quasi-elastic scattering with hyperon";
    case simb::kNCCOH                         : return "neutral current coherent scattering";
    case simb::kCCCOH                         : return "charged current coherent scattering";
    case simb::kNuElectronElastic             : return "electron neutrino elastic";
    case simb::kInverseMuDecay                : return "inverse muon decay";
    default                                   : return "unsupported (" + std::to_string(type) + ")";
  } // switch
} // sim::DumpMCTruth::TruthInteractionTypeName()


//------------------------------------------------------------------------------

std::string sim::DumpMCTruth::ParticleName(int pigid) {
  TParticlePDG const* PDGinfo = TDatabasePDG::Instance()->GetParticle(pigid);
  return PDGinfo? PDGinfo->GetTitle(): ("PDG ID " + std::to_string(pigid));
} // sim::DumpMCTruth::ParticleName()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpMCTruth)

//------------------------------------------------------------------------------
