////////////////////////////////////////////////////////////////////////
// Class:       RecoCheckAna
// Module Type: analyzer
// File:        RecoCheckAna.h
//
// Generated at Fri Jul 15 09:54:26 2011 by Brian Rebel using artmod
// from art v0_07_04.
////////////////////////////////////////////////////////////////////////

#include <map>
#include <set>
#include <utility> // std::move(), std::pair
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Event.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace cheat {
  class RecoCheckAna;
}

class cheat::RecoCheckAna : public art::EDAnalyzer {
public:
  explicit RecoCheckAna(fhicl::ParameterSet const& p);

private:
  void analyze(art::Event const& e) override;
  void beginRun(art::Run const& r) override;

  void CheckReco(
    detinfo::DetectorClocksData const& clockData,
    int const& colID,
    std::vector<recob::Hit const*> const& allhits,
    std::vector<recob::Hit const*> const& colHits,
    std::map<std::pair<int, int>, std::pair<double, double>>& g4RecoBaseIDToPurityEfficiency);
  void CheckRecoClusters(art::Event const& evt,
                         std::string const& label,
                         art::Handle<std::vector<recob::Cluster>> const& clscol,
                         std::vector<recob::Hit const*> const& allhits);
  void CheckRecoTracks(art::Event const& evt,
                       std::string const& label,
                       art::Handle<std::vector<recob::Track>> const& tcol,
                       std::vector<recob::Hit const*> const& allhits);
  void CheckRecoShowers(art::Event const& evt,
                        std::string const& label,
                        art::Handle<std::vector<recob::Shower>> const& scol,
                        std::vector<recob::Hit const*> const& allhits);
  void CheckRecoVertices(art::Event const& evt,
                         std::string const& label,
                         art::Handle<std::vector<recob::Vertex>> const& vtxcol,
                         std::vector<recob::Hit const*> const& allhits);
  void CheckRecoEvents(art::Event const& evt,
                       std::string const& label,
                       art::Handle<std::vector<recob::Event>> const& evtcol,
                       std::vector<recob::Hit const*> const& allhits);
  // method to fill the histograms and TTree
  void FillResults(detinfo::DetectorClocksData const& clockData,
                   std::vector<recob::Hit const*> const& allhits);

  // helper method to the above for clusters, showers and tracks
  void FlattenMap(
    std::map<std::pair<int, int>, std::pair<double, double>> const& g4RecoBaseIDToPurityEfficiency,
    std::map<int, std::vector<std::pair<int, std::pair<double, double>>>>&
      g4IDToRecoBasePurityEfficiency,
    TH1D* purity,
    TH1D* efficiency,
    TH1D* purityEfficiency,
    TH2D* purityEfficiency2D);

  art::ServiceHandle<cheat::BackTrackerService const> fBT;       ///< the back tracker service
  art::ServiceHandle<cheat::ParticleInventoryService const> fPI; ///< the back tracker service

  std::string fHitModuleLabel;     ///< label for module making the hits
  std::string fClusterModuleLabel; ///< label for module making the clusters
  std::string fShowerModuleLabel;  ///< label for module making the showers
  std::string fTrackModuleLabel;   ///< label for module making the tracks
  std::string fVertexModuleLabel;  ///< label for module making the vertices
  std::string fEventModuleLabel;   ///< label for module making the events

  bool fCheckClusters; ///< should we check the reconstruction of clusters?
  bool fCheckShowers;  ///< should we check the reconstruction of showers?
  bool fCheckTracks;   ///< should we check the reconstruction of tracks?
  bool fCheckVertices; ///< should we check the reconstruction of vertices?
  bool fCheckEvents;   ///< should we check the reconstruction of events?

  TH1D* fClusterPurity;             ///< histogram of cluster purity
  TH1D* fClusterEfficiency;         ///< histogram of cluster efficiency
  TH1D* fClusterPurityEfficiency;   ///< histogram of cluster efficiency times purity
  TH2D* fClusterPurityEfficiency2D; ///< scatter histogram of cluster purity and efficiency
  TH1D* fShowerPurity;              ///< histogram of shower purity
  TH1D* fShowerEfficiency;          ///< histogram of shower efficiency
  TH1D* fShowerPurityEfficiency;    ///< histogram of shower efficiency times purity
  TH2D* fShowerPurityEfficiency2D;  ///< scatter histogram of cluster purity and efficiency
  TH1D* fTrackPurity;               ///< histogram of track purity
  TH1D* fTrackEfficiency;           ///< histogram of track efficiency
  TH1D* fTrackPurityEfficiency;     ///< histogram of track efficiency times purity
  TH2D* fTrackPurityEfficiency2D;   ///< scatter histogram of cluster purity and efficiency
  TH1D* fVertexPurity;              ///< histogram of vertex purity
  TH1D* fVertexEfficiency;          ///< histogram of vertex efficiency
  TH1D* fVertexPurityEfficiency;    ///< histogram of vertex efficiency times purity
  TH1D* fEventPurity;               ///< histogram of event purity
  TH1D* fEventEfficiency;           ///< histogram of event efficiency
  TH1D* fEventPurityEfficiency;     ///< histogram of event efficiency times purity

  // The following maps have a pair of the G4 track id and RecoBase object
  // id as the key and then the purity and efficiency (in that order) of the RecoBase object
  // as the value
  std::map<std::pair<int, int>, std::pair<double, double>> fG4ClusterIDToPurityEfficiency;
  std::map<std::pair<int, int>, std::pair<double, double>> fG4ShowerIDToPurityEfficiency;
  std::map<std::pair<int, int>, std::pair<double, double>> fG4TrackIDToPurityEfficiency;

  TTree* fTree;                ///< TTree to save efficiencies
  int frun;                    ///< run number
  int fevent;                  ///< event number
  int ftrackid;                ///< geant track ID
  int fpdg;                    ///< particle pdg code
  double fpmom;                ///< particle momentum
  double fhiteff;              ///< hitfinder efficiency for this particle
  int fnclu;                   ///< number of clusters for this particle
  std::vector<double> fclueff; ///< cluster efficiencies
  std::vector<double> fclupur; ///< cluster purities
  std::vector<int> fcluid;     ///< cluster IDs
  int fnshw;                   ///< number of showers for this particle
  std::vector<double> fshweff; ///< shower efficiencies
  std::vector<double> fshwpur; ///< shower purities
  std::vector<int> fshwid;     ///< shower IDs
  int fntrk;                   ///< number of tracks for this particle
  std::vector<double> ftrkeff; ///< track efficiencies
  std::vector<double> ftrkpur; ///< track purities
  std::vector<int> ftrkid;     ///< track IDs
};

//-------------------------------------------------------------------
cheat::RecoCheckAna::RecoCheckAna(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
  , fHitModuleLabel{p.get<std::string>("HitModuleLabel")}
  , fClusterModuleLabel{p.get<std::string>("ClusterModuleLabel")}
  , fShowerModuleLabel{p.get<std::string>("ShowerModuleLabel")}
  , fTrackModuleLabel{p.get<std::string>("TrackModuleLabel")}
  , fVertexModuleLabel{p.get<std::string>("VertexModuleLabel")}
  , fEventModuleLabel{p.get<std::string>("EventModuleLabel")}
  , fCheckClusters{p.get<bool>("CheckClusters")}
  , fCheckShowers{p.get<bool>("CheckShowers")}
  , fCheckTracks{p.get<bool>("CheckTracks")}
  , fCheckVertices{p.get<bool>("CheckVertices")}
  , fCheckEvents{p.get<bool>("CheckEvents")}
{}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::analyze(art::Event const& e)
{
  // check that this is MC, stop if it isn't
  if (e.isRealData()) {
    mf::LogWarning("RecoVetter") << "attempting to run MC truth check on "
                                 << "real data, bail";
    return;
  }

  // get all hits in the event to figure out how many there are
  auto const hithdl = e.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);
  std::vector<recob::Hit const*> allhits;
  allhits.reserve(hithdl->size());
  for (auto const& hit : *hithdl)
    allhits.push_back(&hit);

  if (fCheckClusters) {
    auto const clscol = e.getHandle<std::vector<recob::Cluster>>(fClusterModuleLabel);
    if (clscol) this->CheckRecoClusters(e, fClusterModuleLabel, clscol, allhits);
  }
  if (fCheckTracks) {
    auto const trkcol = e.getHandle<std::vector<recob::Track>>(fTrackModuleLabel);
    if (trkcol) this->CheckRecoTracks(e, fTrackModuleLabel, trkcol, allhits);
  }
  if (fCheckShowers) {
    auto const shwcol = e.getHandle<std::vector<recob::Shower>>(fShowerModuleLabel);
    if (shwcol) this->CheckRecoShowers(e, fShowerModuleLabel, shwcol, allhits);
  }
  if (fCheckVertices) {
    auto const vtxcol = e.getHandle<std::vector<recob::Vertex>>(fVertexModuleLabel);
    if (vtxcol) this->CheckRecoVertices(e, fVertexModuleLabel, vtxcol, allhits);
  }
  if (fCheckEvents) {
    auto const evtcol = e.getHandle<std::vector<recob::Event>>(fEventModuleLabel);
    if (evtcol) this->CheckRecoEvents(e, fEventModuleLabel, evtcol, allhits);
  }

  frun = e.run();
  fevent = e.id().event();

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  this->FillResults(clockData, allhits);

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::beginRun(art::Run const& /*r*/)
{
  art::ServiceHandle<art::TFileService const> tfs;

  if (fCheckEvents) {
    fEventPurity = tfs->make<TH1D>("eventPurity", ";Purity;Events", 100, 0., 1.1);
    fEventEfficiency = tfs->make<TH1D>("eventEfficiency", ";Efficiency;Events", 100, 0., 1.1);
    fEventPurityEfficiency =
      tfs->make<TH1D>("eventPurityEfficiency", ";purityEfficiency;Events", 110, 0., 1.1);
  }
  if (fCheckVertices) {
    fVertexPurity = tfs->make<TH1D>("vertexPurity", ";Purity;Vertices", 100, 0., 1.1);
    fVertexEfficiency = tfs->make<TH1D>("vertexEfficiency", ";Efficiency;Vertices", 100, 0., 1.1);
    fVertexPurityEfficiency =
      tfs->make<TH1D>("vertexPurityEfficiency", ";purityEfficiency;Vertex", 110, 0., 1.1);
  }
  if (fCheckTracks) {
    fTrackPurity = tfs->make<TH1D>("trackPurity", ";Purity;Tracks", 100, 0., 1.1);
    fTrackEfficiency = tfs->make<TH1D>("trackEfficiency", ";Efficiency;Tracks", 100, 0., 1.1);
    fTrackPurityEfficiency =
      tfs->make<TH1D>("trackPurityEfficiency", ";purityEfficiency;Tracks", 110, 0., 1.1);
    fTrackPurityEfficiency2D =
      tfs->make<TH2D>("trackPurityEfficiency2D", ";purity;efficiency", 110, 0., 1.1, 110, 0., 1.1);
  }
  if (fCheckShowers) {
    fShowerPurity = tfs->make<TH1D>("showerPurity", ";Purity;Showers", 100, 0., 1.1);
    fShowerEfficiency = tfs->make<TH1D>("showerEfficiency", ";Efficiency;Showers", 100, 0., 1.1);
    fShowerPurityEfficiency =
      tfs->make<TH1D>("showerPurityEfficiency", ";purityEfficiency;Showers", 110, 0., 1.1);
    fShowerPurityEfficiency2D =
      tfs->make<TH2D>("showerPurityEfficiency2D", ";purity;efficiency", 110, 0., 1.1, 110, 0., 1.1);
  }
  if (fCheckClusters) {
    fClusterPurity = tfs->make<TH1D>("clusterPurity", ";Purity;Clusters", 110, 0., 1.1);
    fClusterEfficiency = tfs->make<TH1D>("clusterEfficiency", ";Efficiency;Clusters", 110, 0., 1.1);
    fClusterPurityEfficiency =
      tfs->make<TH1D>("clusterPurityEfficiency", ";purityEfficiency;Clusters", 110, 0., 1.1);
    fClusterPurityEfficiency2D = tfs->make<TH2D>(
      "clusterPurityEfficiency2D", ";purity;efficiency", 110, 0., 1.1, 110, 0., 1.1);
  }

  fTree = tfs->make<TTree>("cheatertree", "cheater tree");
  fTree->Branch("run", &frun, "run/I");
  fTree->Branch("event", &fevent, "event/I");
  fTree->Branch("trackid", &ftrackid, "trackid/I");
  fTree->Branch("pdg", &fpdg, "pdg/I");
  fTree->Branch("pmom", &fpmom, "pmom/D");
  fTree->Branch("hiteff", &fhiteff, "hiteff/D");
  fTree->Branch("nclu", &fnclu, "nclu/I");
  fTree->Branch("clueff", &fclueff);
  fTree->Branch("clupur", &fclupur);
  fTree->Branch("cluid", &fcluid);
  fTree->Branch("nshw", &fnshw, "nshw/I");
  fTree->Branch("shweff", &fshweff);
  fTree->Branch("shwpur", &fshwpur);
  fTree->Branch("shwid", &fshwid);
  fTree->Branch("ntrk", &fntrk, "ntrk/I");
  fTree->Branch("trkeff", &ftrkeff);
  fTree->Branch("trkpur", &ftrkpur);
  fTree->Branch("trkid", &ftrkid);

  return;
}

//-------------------------------------------------------------------
// colID is the ID of the RecoBase object and colHits are the recob::Hits
// associated with it
void cheat::RecoCheckAna::CheckReco(
  detinfo::DetectorClocksData const& clockData,
  int const& colID,
  std::vector<recob::Hit const*> const& allhits,
  std::vector<recob::Hit const*> const& colHits,
  std::map<std::pair<int, int>, std::pair<double, double>>& g4RecoBaseIDToPurityEfficiency)
{
  auto const view = colHits[0]->View();

  for (int trackID : fBT->GetSetOfTrackIds(clockData, colHits)) {
    double purity = fBT->HitCollectionPurity(clockData, {trackID}, colHits);
    double efficiency = fBT->HitCollectionEfficiency(clockData, {trackID}, colHits, allhits, view);
    g4RecoBaseIDToPurityEfficiency[{trackID, colID}] = {purity, efficiency};
  }

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::CheckRecoClusters(art::Event const& evt,
                                            std::string const& label,
                                            art::Handle<std::vector<recob::Cluster>> const& clscol,
                                            std::vector<recob::Hit const*> const& allhits)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  art::FindMany<recob::Hit> fmh(clscol, evt, label);

  for (size_t c = 0; c < clscol->size(); ++c)
    this->CheckReco(clockData, clscol->at(c).ID(), allhits, fmh.at(c), fG4ClusterIDToPurityEfficiency);

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::CheckRecoTracks(art::Event const& evt,
                                          std::string const& label,
                                          art::Handle<std::vector<recob::Track>> const& tcol,
                                          std::vector<recob::Hit const*> const& allhits)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  art::FindMany<recob::Hit> fmh(tcol, evt, label);

  for (size_t p = 0; p < tcol->size(); ++p)
    this->CheckReco(clockData, tcol->at(p).ID(), allhits, fmh.at(p), fG4TrackIDToPurityEfficiency);

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::CheckRecoShowers(art::Event const& evt,
                                           std::string const& label,
                                           art::Handle<std::vector<recob::Shower>> const& scol,
                                           std::vector<recob::Hit const*> const& allhits)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  art::FindMany<recob::Hit> fmh(scol, evt, label);

  for (size_t p = 0; p < scol->size(); ++p)
    this->CheckReco(clockData, scol->at(p).ID(), allhits, fmh.at(p), fG4ShowerIDToPurityEfficiency);

  return;
}

//-------------------------------------------------------------------
//a true vertex will either consist of primary particles originating from
//the interaction vertex, or a primary particle decaying to make daughters
void cheat::RecoCheckAna::CheckRecoVertices(art::Event const& evt,
                                            std::string const& label,
                                            art::Handle<std::vector<recob::Vertex>> const& vtxcol,
                                            std::vector<recob::Hit const*> const& allhits)
{
  const sim::ParticleList& plist = fPI->ParticleList();

  std::vector<std::set<int>> ids(1);
  // loop over all primary particles and put their ids into the first set of the
  // vector.  add another set for each primary particle that also has daughters
  // and put those daughters into the new set
  // PartPair is a (track ID, particle pointer) pair
  for (const auto& PartPair : plist) {
    auto trackID = PartPair.first;
    if (!plist.IsPrimary(trackID)) continue;
    const simb::MCParticle& part = *(PartPair.second);
    ids[0].insert(trackID);
    if (part.NumberDaughters() > 0) {
      std::set<int> dv;
      for (int d = 0; d < part.NumberDaughters(); ++d)
        dv.insert(part.Daughter(d));
      ids.push_back(std::move(dv));
    } // end if this primary particle has daughters
  }   // end loop over primaries

  art::FindMany<recob::Hit> fmh(vtxcol, evt, label);

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  for (size_t v = 0; v < vtxcol->size(); ++v) {
    auto const& hits = fmh.at(v);

    double maxPurity = -1.;
    double maxEfficiency = -1.;

    for (auto const& idSet : ids) {
      double purity = fBT->HitCollectionPurity(clockData, idSet, hits);
      double efficiency = fBT->HitCollectionEfficiency(clockData, idSet, hits, allhits, geo::k3D);

      if (purity > maxPurity) maxPurity = purity;
      if (efficiency > maxEfficiency) maxEfficiency = efficiency;
    }

    fVertexPurity->Fill(maxPurity);
    fVertexEfficiency->Fill(maxEfficiency);
    fVertexPurityEfficiency->Fill(maxPurity * maxEfficiency);

  } // end loop over vertices

  return;
}

//-------------------------------------------------------------------
// in this method one should loop over the primary particles from a given
// MCTruth collection
/// \todo need to divy it up in the case where there is more than 1 true interaction in a spill
void cheat::RecoCheckAna::CheckRecoEvents(art::Event const& evt,
                                          std::string const& label,
                                          art::Handle<std::vector<recob::Event>> const& evtcol,
                                          std::vector<recob::Hit const*> const& allhits)
{
  const sim::ParticleList& plist = fPI->ParticleList();

  // loop over all primaries in the plist and grab them and their daughters to put into
  // the set of track ids to pass on to the back tracker
  std::set<int> ids;
  for (const auto& PartPair : plist) {
    auto trackID = PartPair.first;
    if (!plist.IsPrimary(trackID)) continue;
    const simb::MCParticle& part = *(PartPair.second);
    ids.insert(trackID);
    for (int d = 0; d < part.NumberDaughters(); ++d)
      ids.insert(part.Daughter(d));
  } // end loop over primaries

  art::FindMany<recob::Hit> fmh(evtcol, evt, label);

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  for (size_t ev = 0; ev < evtcol->size(); ++ev) {
    auto const& hits = fmh.at(ev);

    double purity = fBT->HitCollectionPurity(clockData, ids, hits);
    double efficiency = fBT->HitCollectionEfficiency(clockData, ids, hits, allhits, geo::k3D);

    fEventPurity->Fill(purity);
    fEventEfficiency->Fill(efficiency);
    fEventPurityEfficiency->Fill(purity * efficiency);

  } // end loop over events

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::FlattenMap(
  std::map<std::pair<int, int>, std::pair<double, double>> const& g4RecoBaseIDToPurityEfficiency,
  std::map<int, std::vector<std::pair<int, std::pair<double, double>>>>&
    g4IDToRecoBasePurityEfficiency,
  TH1D* purity,
  TH1D* efficiency,
  TH1D* purityEfficiency,
  TH2D* purityEfficiency2D)
{

  // map of key cluster ID to pair of purity, efficiency
  std::map<int, std::pair<double, double>> recoBIDToPurityEfficiency;

  for (auto const& [g4cl, pe] : g4RecoBaseIDToPurityEfficiency) {
    auto const& [trackID, recoID] = g4cl;

    // g4IDToRecoBasePurityEfficiency is a map with key of trackID of a vector of clusterIDs of pairs of purity and efficiency
    g4IDToRecoBasePurityEfficiency[trackID].emplace_back(recoID, pe);

    // now find the maximum purity to determine the purity and efficiency
    // for this RecoBase object
    auto it = recoBIDToPurityEfficiency.find(recoID);
    if (it != recoBIDToPurityEfficiency.end()) {
      if (pe.first > it->second.first) it->second = pe;
    }
    else
      recoBIDToPurityEfficiency[recoID] = pe;
  }

  // now fill the histograms
  for (auto const& [recoID, pe] : recoBIDToPurityEfficiency) {
    purity->Fill(pe.first);
    efficiency->Fill(pe.second);
    purityEfficiency->Fill(pe.first * pe.second);
    purityEfficiency2D->Fill(pe.first, pe.second);
  }

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::FillResults(detinfo::DetectorClocksData const& clockData,
                                      std::vector<recob::Hit const*> const& allhits)
{
  // map the g4 track id to energy deposited in a hit
  std::map<int, double> g4IDToHitEnergy;
  for (auto const* hit : allhits) {
    for (auto const& ide : fBT->HitToTrackIDEs(clockData, *hit))
      g4IDToHitEnergy[ide.trackID] += ide.energy;
  } // end loop over hits to fill map

  // flatten the G4RecoBaseIDToPurityEfficiency maps to have just the g4ID as
  // the key and the rest of the information in vector form
  std::map<int, std::vector<std::pair<int, std::pair<double, double>>>>
    g4IDToClusterPurityEfficiency;
  std::map<int, std::vector<std::pair<int, std::pair<double, double>>>>
    g4IDToShowerPurityEfficiency;
  std::map<int, std::vector<std::pair<int, std::pair<double, double>>>> g4IDToTrackPurityEfficiency;

  if (fCheckClusters)
    this->FlattenMap(fG4ClusterIDToPurityEfficiency,
                     g4IDToClusterPurityEfficiency,
                     fClusterPurity,
                     fClusterEfficiency,
                     fClusterPurityEfficiency,
                     fClusterPurityEfficiency2D);
  if (fCheckShowers)
    this->FlattenMap(fG4ShowerIDToPurityEfficiency,
                     g4IDToShowerPurityEfficiency,
                     fShowerPurity,
                     fShowerEfficiency,
                     fShowerPurityEfficiency,
                     fShowerPurityEfficiency2D);
  if (fCheckTracks)
    this->FlattenMap(fG4TrackIDToPurityEfficiency,
                     g4IDToTrackPurityEfficiency,
                     fTrackPurity,
                     fTrackEfficiency,
                     fTrackPurityEfficiency,
                     fTrackPurityEfficiency2D);

  // fill the tree vectors
  for (int trackID : fBT->GetSetOfTrackIds()) {

    const simb::MCParticle* part = fPI->TrackIdToParticle(trackID);

    ftrackid = std::abs(trackID);
    fpdg = part->PdgCode();
    fpmom = part->P();

    // figure out how much of the energy deposited from this particle is stored in hits
    double totalDep = 0.;
    for (auto const* ide : fBT->TrackIdToSimIDEs_Ps(trackID))
      totalDep += ide->energy;

    if (totalDep > 0.) fhiteff = g4IDToHitEnergy[trackID] / totalDep;

    std::vector<std::pair<int, std::pair<double, double>>> clVec;
    std::vector<std::pair<int, std::pair<double, double>>> shVec;
    std::vector<std::pair<int, std::pair<double, double>>> trVec;

    if (auto it = g4IDToClusterPurityEfficiency.find(trackID);
        it != g4IDToClusterPurityEfficiency.end())
      clVec = it->second;

    if (auto it = g4IDToShowerPurityEfficiency.find(trackID);
        it != g4IDToShowerPurityEfficiency.end())
      shVec = it->second;

    if (auto it = g4IDToTrackPurityEfficiency.find(trackID);
        it != g4IDToTrackPurityEfficiency.end())
      trVec = it->second;

    fnclu = clVec.size();
    fnshw = shVec.size();
    fntrk = trVec.size();

    for (auto const& [id, pe] : clVec) {
      fcluid.push_back(id);
      fclupur.push_back(pe.first);
      fclueff.push_back(pe.second);
    }

    for (auto const& [id, pe] : shVec) {
      fshwid.push_back(id);
      fshwpur.push_back(pe.first);
      fshweff.push_back(pe.second);
    }

    for (auto const& [id, pe] : trVec) {
      ftrkid.push_back(id);
      ftrkpur.push_back(pe.first);
      ftrkeff.push_back(pe.second);
    }

    fTree->Fill();
  }

  // clean up for the next event

  // clear the maps of G4 track id to efficiency and purity for
  // various RecoBase objects
  fG4ClusterIDToPurityEfficiency.clear();
  fG4ShowerIDToPurityEfficiency.clear();
  fG4TrackIDToPurityEfficiency.clear();

  // clear the vectors hooked up to the tree
  fclueff.clear();
  fclupur.clear();
  fcluid.clear();
  ftrkeff.clear();
  ftrkpur.clear();
  ftrkid.clear();
  fshweff.clear();
  fshwpur.clear();
  fshwid.clear();

  return;
}

DEFINE_ART_MODULE(cheat::RecoCheckAna)
