/**
 * @file   DumpMCParticles_module.cc
 * @brief  Module dumping MCarticle information on screen
 * @date   December 3rd, 2015
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 */

// LArSoft libraries
#include "lardataalg/MCDumpers/MCDumpers.h" // sim::dump namespace
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

// nusimdata libraries
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <string>

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

    fhicl::OptionalAtom<art::InputTag> ParticleTruthInfo {
      Name("ParticleTruthInfo"),
      Comment
        ("label of the association to MCTruth (default: as `InputParticles`)")
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

  /// May print some warnings.
  void endJob() override;

  /**
   * @brief Dumps the content of the specified particle in the output stream.
   * @tparam Stream the type of output stream
   * @param out the output stream
   * @param particle the particle to be dumped
   * @param truthTag input tag of the truth record the particle derived from
   * @param truthIndex index of particle in the truth record this derived from
   * @param indent base indentation string (default: none)
   * @param bIndentFirst if first output line should be indented (default: yes)
   *
   * The indent string is prepended to every line of output, with the possible
   * exception of the first one, in case bIndentFirst is true.
   *
   * If `truthTag` module label is empty, it is assumed that this information
   * could not be retrieved, and it will be silently omitted.
   * If `truthIndex` is `sim::NoGeneratorIndex`, it is assumed that this
   * information could not be retrieved, and it will be silently omitted.
   *
   * The output starts on the current line, and the last line is NOT broken.
   */
  template <typename Stream>
  void DumpMCParticle(
    Stream&& out, simb::MCParticle const& particle,
    art::InputTag const& truthTag, sim::GeneratedParticleInfo const& truthInfo,
    std::string indent = "", bool bIndentFirst = true
    ) const;


    private:

  art::InputTag fInputParticles; ///< name of MCParticle's data product
  art::InputTag fParticleTruthInfo; ///< name of MCParticle assns data product
  std::string fOutputCategory; ///< name of the stream for output
  unsigned int fPointsPerLine; ///< trajectory points per output line

  unsigned int fNEvents = 0U; ///< Count of processed events.
  /// Count of events without truth association.
  unsigned int fNMissingTruth = 0U;
  /// Count of events without truth index.
  unsigned int fNMissingTruthIndex = 0U;

}; // class sim::DumpMCParticles


//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
namespace {

  //----------------------------------------------------------------------------
  class ProductNameCache {

      public:
    ProductNameCache(art::Event const& event): fEvent(event) {}

    /// Returns the input tag corresponding to the specified art product id.
    template <typename T>
    art::InputTag const& operator[](art::Ptr<T> const& ptr)
      {
        auto const iInfo = fNames.find(ptr.id());
        return (iInfo == fNames.end())? fetch(ptr): iInfo->second;
      }

      private:
    art::Event const& fEvent;
    std::map<art::ProductID, art::InputTag> fNames;

    template <typename T>
    art::InputTag fetchTag(art::Ptr<T> const& ptr)
      {
        art::Handle<std::vector<T>> handle;
        return fEvent.get(ptr.id(), handle)
          ? handle.provenance()->inputTag()
          : art::InputTag{}
          ;
      }

    template <typename T>
    art::InputTag const& fetch(art::Ptr<T> const& ptr)
      {
        art::InputTag const tag = fetchTag(ptr);
        return fNames.emplace(ptr.id(), tag).first->second;
      }

  }; // class ProductNameCache


  //----------------------------------------------------------------------------
  template <typename Right, typename Metadata, typename Left>
  std::unique_ptr<art::FindOneP<Right, Metadata>> makeFindOneP(
    art::ValidHandle<std::vector<Left>> const& handle,
    art::Event const& event,
    art::InputTag const& tag
  ) {

    art::Handle<art::Assns<Left, Right, Metadata>> assnsHandle;
    if (!event.getByLabel(tag, assnsHandle)) return {};

    return std::make_unique<art::FindOneP<Right, Metadata>>
      (handle, event, tag);

  } // makeFindOneP()


  //----------------------------------------------------------------------------

} // local namespace


//------------------------------------------------------------------------------
sim::DumpMCParticles::DumpMCParticles(Parameters const& config)
  : EDAnalyzer(config)
  , fInputParticles(config().InputParticles())
  , fOutputCategory(config().OutputCategory())
  , fPointsPerLine(config().PointsPerLine())
{
  if (!config().ParticleTruthInfo(fParticleTruthInfo))
    fParticleTruthInfo = fInputParticles;

}

//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpMCParticles::DumpMCParticle(
  Stream&& out, simb::MCParticle const& particle,
  art::InputTag const& truthTag, sim::GeneratedParticleInfo const& truthInfo,
  std::string indent /* = "" */, bool bIndentFirst /* = true */
) const {

  if (!truthTag.label().empty() || truthInfo.hasGeneratedParticleIndex()) {
    out << "(from ";
    if (truthTag.label().empty()) out << "unknown truth record";
    else out << "'" << truthTag.encode() << "'";
    if (truthInfo.hasGeneratedParticleIndex())
      out << " particle #" << truthInfo.generatedParticleIndex();
    out << ") ";
  }

  sim::dump::DumpMCParticle
    (std::forward<Stream>(out), particle, indent, bIndentFirst? indent: "");

  const unsigned int nPoints = particle.NumberTrajectoryPoints();
  if ((nPoints > 0) && (fPointsPerLine > 0)) {
    out << ":";
    sim::dump::DumpMCParticleTrajectory(
      std::forward<Stream>(out), particle.Trajectory(),
      fPointsPerLine, indent + "  "
      );
  } // if has points

} // sim::DumpMCParticles::DumpMCParticle()


//------------------------------------------------------------------------------
void sim::DumpMCParticles::analyze(art::Event const& event) {

  ++fNEvents;

  ProductNameCache namesRegistry(event);

  // get the particles from the event
  auto const& particleHandle
    = event.getValidHandle<std::vector<simb::MCParticle>>(fInputParticles);
  auto const& Particles = *particleHandle;

  // get the association to MCTruth
  // - try first the more complete one, with true particle indices
  // - as a fallback, go without true particle indices
  auto particleToTruth = makeFindOneP<simb::MCTruth, sim::GeneratedParticleInfo>
    (particleHandle, event, fParticleTruthInfo);
  std::unique_ptr<art::FindOneP<simb::MCTruth>> particleToTruthLight;
  if (!particleToTruth) {
    ++fNMissingTruthIndex;
    particleToTruthLight = makeFindOneP<simb::MCTruth, void>
      (particleHandle, event, fParticleTruthInfo);
    if (!particleToTruthLight) ++fNMissingTruth;
  }

  mf::LogVerbatim(fOutputCategory) << "Event " << event.id()
    << ": data product '" << fInputParticles.encode() << "' contains "
    << Particles.size() << " MCParticle's";

  unsigned int iParticle = 0;
  for (simb::MCParticle const& particle: Particles) {
    // flush on every particle,
    // since the output buffer might grow too large otherwise
    mf::LogVerbatim log(fOutputCategory);

    // fetch the input tag of the truth information (if any)
    art::Ptr<simb::MCTruth> const& truth = particleToTruth
      ? particleToTruth->at(iParticle)
      : particleToTruthLight
        ? particleToTruthLight->at(iParticle)
        : art::Ptr<simb::MCTruth>{}
      ;
    art::InputTag const& truthTag
      = truth? namesRegistry[truth]: art::InputTag{};

    // fetch the index of the true particle in the truth record (if any)
    sim::GeneratedParticleInfo truthInfo = particleToTruth
      ? particleToTruth->data(iParticle).ref()
      : sim::GeneratedParticleInfo::NoGeneratedParticleIndex
      ;

    // a bit of a header
    log << "\n[#" << (iParticle++) << "] ";
    DumpMCParticle(log, particle, truthTag, truthInfo, "  ", false);
  } // for

  mf::LogVerbatim(fOutputCategory) << "\n";

} // sim::DumpMCParticles::analyze()


//------------------------------------------------------------------------------
void sim::DumpMCParticles::endJob() {

  if (fNMissingTruth > 0) {
    mf::LogProblem(fOutputCategory)
      << "Warning: " << fNMissingTruth << "/" << fNEvents
      << " events lacked event generator information for '"
      << fParticleTruthInfo << "'.";
  }
  else if (fNMissingTruthIndex > 0) {
    mf::LogProblem(fOutputCategory)
      << "Warning: " << fNMissingTruthIndex << "/" << fNEvents
      << " events lacked information of which particles of '"
      << fParticleTruthInfo << "' are generator particles.";
  }

} // sim::DumpMCParticles::endJob()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpMCParticles)

//------------------------------------------------------------------------------
