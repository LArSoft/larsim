////////////////////////////////////////////////////////////////////////
// Class:       RecoCheckAna
// Module Type: analyzer
// File:        RecoCheckAna.h
//
// Generated at Fri Jul 15 09:54:26 2011 by Brian Rebel using artmod
// from art v0_07_04.
////////////////////////////////////////////////////////////////////////

#include <utility> // std::move()

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Event.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nug4/ParticleNavigation/ParticleList.h"

namespace cheat {
  class RecoCheckAna;
}

class TH1D;
class TTree;

class cheat::RecoCheckAna : public art::EDAnalyzer {
public:
  explicit RecoCheckAna(fhicl::ParameterSet const &p);

  virtual void analyze(art::Event const &e);

  virtual void beginRun(art::Run const &r);



private:

  void CheckReco        (int                                                 const& colID,
			 std::vector< art::Ptr<recob::Hit> > 	 	     const& allhits,
			 std::vector< art::Ptr<recob::Hit> > 	 	     const& colHits,
			 std::map<std::pair<int, int>, std::pair<double, double> >& g4RecoBaseIDToPurityEfficiency);
  void CheckRecoClusters(art::Event                                 const& evt,
			 std::string                                const& label,
			 art::Handle< std::vector<recob::Cluster> > const& clscol,
			 std::vector< art::Ptr<recob::Hit> >        const& allhits);
  void CheckRecoTracks  (art::Event                                 const& evt,
			 std::string                                const& label,
			 art::Handle< std::vector<recob::Track> >   const& tcol,
			 std::vector< art::Ptr<recob::Hit> >        const& allhits);
  void CheckRecoShowers (art::Event                                const& evt,
			 std::string                               const& label,
			 art::Handle< std::vector<recob::Shower> > const& scol,
			 std::vector< art::Ptr<recob::Hit> >       const& allhits);
  void CheckRecoVertices(art::Event                                 const& evt,
			 std::string                                const& label,
			 art::Handle< std::vector<recob::Vertex> >  const& vtxcol,
			 std::vector< art::Ptr<recob::Hit> >        const& allhits);
  void CheckRecoEvents  (art::Event                                 const& evt,
		         std::string                                const& label,
		         art::Handle< std::vector<recob::Event> >   const& evtcol,
		         std::vector< art::Ptr<recob::Hit> >        const& allhits);
  // method to fill the histograms and TTree
  void FillResults(std::vector< art::Ptr<recob::Hit> > const& allhits);

  // helper method to the above for clusters, showers and tracks
  void FlattenMap(std::map<std::pair<int, int>,  std::pair<double, double> > const& g4RecoBaseIDToPurityEfficiency,
		  std::map<int, std::vector<std::pair<int, std::pair<double, double> > > >& g4IDToRecoBasePurityEfficiency,
                  TH1D*                                                            purity,
		  TH1D*                                                            efficiency,
		  TH1D*                                                            purityEfficiency,
		  TH2D*                                                            purityEfficiency2D);

  art::ServiceHandle<cheat::BackTrackerService const> fBT; ///< the back tracker service
  art::ServiceHandle<cheat::ParticleInventoryService const> fPI; ///< the back tracker service

  std::string fHitModuleLabel;		 ///< label for module making the hits
  std::string fClusterModuleLabel;	 ///< label for module making the clusters
  std::string fShowerModuleLabel; 	 ///< label for module making the showers
  std::string fTrackModuleLabel;  	 ///< label for module making the tracks
  std::string fVertexModuleLabel; 	 ///< label for module making the vertices
  std::string fEventModuleLabel;  	 ///< label for module making the events

  bool        fCheckClusters;     	 ///< should we check the reconstruction of clusters?
  bool        fCheckShowers;      	 ///< should we check the reconstruction of showers?
  bool        fCheckTracks;       	 ///< should we check the reconstruction of tracks?
  bool        fCheckVertices;     	 ///< should we check the reconstruction of vertices?
  bool        fCheckEvents;       	 ///< should we check the reconstruction of events?

  TH1D*       fClusterPurity;     	 ///< histogram of cluster purity
  TH1D*       fClusterEfficiency; 	 ///< histogram of cluster efficiency
  TH1D*       fClusterPurityEfficiency;  ///< histogram of cluster efficiency times purity
  TH2D*       fClusterPurityEfficiency2D;///< scatter histogram of cluster purity and efficiency
  TH1D*       fShowerPurity;		 ///< histogram of shower purity
  TH1D*       fShowerEfficiency;  	 ///< histogram of shower efficiency
  TH1D*       fShowerPurityEfficiency;   ///< histogram of shower efficiency times purity
  TH2D*       fShowerPurityEfficiency2D; ///< scatter histogram of cluster purity and efficiency
  TH1D*       fTrackPurity;       	 ///< histogram of track purity
  TH1D*       fTrackEfficiency;   	 ///< histogram of track efficiency
  TH1D*       fTrackPurityEfficiency;    ///< histogram of track efficiency times purity
  TH2D*       fTrackPurityEfficiency2D;  ///< scatter histogram of cluster purity and efficiency
  TH1D*       fVertexPurity;      	 ///< histogram of vertex purity
  TH1D*       fVertexEfficiency;  	 ///< histogram of vertex efficiency
  TH1D*       fVertexPurityEfficiency;   ///< histogram of vertex efficiency times purity
  TH1D*       fEventPurity;       	 ///< histogram of event purity
  TH1D*       fEventEfficiency;   	 ///< histogram of event efficiency
  TH1D*       fEventPurityEfficiency;    ///< histogram of event efficiency times purity

  // The following maps have a pair of the G4 track id and RecoBase object
  // id as the key and then the purity and efficiency (in that order) of the RecoBase object
  // as the value
  std::map< std::pair<int, int>, std::pair<double, double> > fG4ClusterIDToPurityEfficiency;
  std::map< std::pair<int, int>, std::pair<double, double> > fG4ShowerIDToPurityEfficiency;
  std::map< std::pair<int, int>, std::pair<double, double> > fG4TrackIDToPurityEfficiency;

  TTree*              fTree;        ///< TTree to save efficiencies
  int                 frun;         ///< run number
  int                 fevent;       ///< event number
  int                 ftrackid;     ///< geant track ID
  int                 fpdg;         ///< particle pdg code
  double              fpmom;        ///< particle momentum
  double              fhiteff;      ///< hitfinder efficiency for this particle
  int                 fnclu;        ///< number of clusters for this particle
  std::vector<double> fclueff;      ///< cluster efficiencies
  std::vector<double> fclupur;      ///< cluster purities
  std::vector<int>    fcluid;       ///< cluster IDs
  int                 fnshw;        ///< number of showers for this particle
  std::vector<double> fshweff;      ///< shower efficiencies
  std::vector<double> fshwpur;      ///< shower purities
  std::vector<int>    fshwid;       ///< shower IDs
  int                 fntrk;        ///< number of tracks for this particle
  std::vector<double> ftrkeff;      ///< track efficiencies
  std::vector<double> ftrkpur;      ///< track purities
  std::vector<int>    ftrkid;       ///< track IDs

};

//-------------------------------------------------------------------
cheat::RecoCheckAna::RecoCheckAna(fhicl::ParameterSet const &p)
  : EDAnalyzer(p)
  , fHitModuleLabel    {p.get< std::string >("HitModuleLabel")}
  , fClusterModuleLabel{p.get< std::string >("ClusterModuleLabel")}
  , fShowerModuleLabel {p.get< std::string >("ShowerModuleLabel"  )}
  , fTrackModuleLabel  {p.get< std::string >("TrackModuleLabel" )}
  , fVertexModuleLabel {p.get< std::string >("VertexModuleLabel" )}
  , fEventModuleLabel  {p.get< std::string >("EventModuleLabel"  )}
  , fCheckClusters     {p.get< bool        >("CheckClusters")}
  , fCheckShowers      {p.get< bool        >("CheckShowers" )}
  , fCheckTracks       {p.get< bool        >("CheckTracks"  )}
  , fCheckVertices     {p.get< bool        >("CheckVertices")}
  , fCheckEvents       {p.get< bool        >("CheckEvents"  )}
{}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::analyze(art::Event const &e)
{
  // check that this is MC, stop if it isn't
  if(e.isRealData()){
    mf::LogWarning("RecoVetter") << "attempting to run MC truth check on "
				 << "real data, bail";
    return;
  }

  // get all hits in the event to figure out how many there are
  art::Handle< std::vector<recob::Hit> > hithdl;
  e.getByLabel(fHitModuleLabel, hithdl);
  std::vector< art::Ptr<recob::Hit> > allhits;
  art::fill_ptr_vector(allhits, hithdl);

  // define variables to hold the reconstructed objects
  art::Handle< std::vector<recob::Cluster> > clscol;
  art::Handle< std::vector<recob::Track> >   trkcol;
  art::Handle< std::vector<recob::Shower> >  shwcol;
  art::Handle< std::vector<recob::Vertex> >  vtxcol;
  art::Handle< std::vector<recob::Event> >   evtcol;

  if(fCheckClusters){
    e.getByLabel(fClusterModuleLabel, clscol);
    if( !clscol.failedToGet() ) this->CheckRecoClusters(e, fClusterModuleLabel, clscol, allhits);
  }
  if(fCheckTracks){
    e.getByLabel(fTrackModuleLabel, trkcol);
    if( !trkcol.failedToGet() ) this->CheckRecoTracks(e, fTrackModuleLabel, trkcol, allhits);
  }
  if(fCheckShowers){
    e.getByLabel(fShowerModuleLabel, shwcol);
    if( !shwcol.failedToGet() ) this->CheckRecoShowers(e, fShowerModuleLabel, shwcol, allhits);
  }
  if(fCheckVertices){
    e.getByLabel(fVertexModuleLabel, vtxcol);
    if( !vtxcol.failedToGet() ) this->CheckRecoVertices(e, fVertexModuleLabel, vtxcol, allhits);
  }
  if(fCheckEvents){
    e.getByLabel(fEventModuleLabel, evtcol);
    if( !evtcol.failedToGet() ) this->CheckRecoEvents(e, fEventModuleLabel, evtcol, allhits);
  }

  frun = e.run();
  fevent = e.id().event();

  this->FillResults(allhits);

  return;

}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::beginRun(art::Run const &/*r*/)
{
  art::ServiceHandle<art::TFileService const> tfs;

  if(fCheckEvents){
    fEventPurity	     = tfs->make<TH1D>("eventPurity",       ";Purity;Events",       100, 0., 1.1);
    fEventEfficiency  	     = tfs->make<TH1D>("eventEfficiency",   ";Efficiency;Events",   100, 0., 1.1);
    fEventPurityEfficiency   = tfs->make<TH1D>("eventPurityEfficiency", ";purityEfficiency;Events", 110, 0., 1.1);
  }
  if(fCheckVertices){
    fVertexPurity     	     = tfs->make<TH1D>("vertexPurity",      ";Purity;Vertices",     100, 0., 1.1);
    fVertexEfficiency 	     = tfs->make<TH1D>("vertexEfficiency",  ";Efficiency;Vertices", 100, 0., 1.1);
    fVertexPurityEfficiency  = tfs->make<TH1D>("vertexPurityEfficiency", ";purityEfficiency;Vertex", 110, 0., 1.1);
  }
  if(fCheckTracks){
    fTrackPurity      	     = tfs->make<TH1D>("trackPurity",       ";Purity;Tracks",       100, 0., 1.1);
    fTrackEfficiency  	     = tfs->make<TH1D>("trackEfficiency",   ";Efficiency;Tracks",   100, 0., 1.1);
    fTrackPurityEfficiency   = tfs->make<TH1D>("trackPurityEfficiency", ";purityEfficiency;Tracks", 110, 0., 1.1);
    fTrackPurityEfficiency2D = tfs->make<TH2D>("trackPurityEfficiency2D", ";purity;efficiency", 110, 0., 1.1, 110, 0., 1.1);
  }
  if(fCheckShowers){
    fShowerPurity     	     = tfs->make<TH1D>("showerPurity",      ";Purity;Showers",      100, 0., 1.1);
    fShowerEfficiency 	     = tfs->make<TH1D>("showerEfficiency",  ";Efficiency;Showers",  100, 0., 1.1);
    fShowerPurityEfficiency  = tfs->make<TH1D>("showerPurityEfficiency", ";purityEfficiency;Showers", 110, 0., 1.1);
    fShowerPurityEfficiency2D= tfs->make<TH2D>("showerPurityEfficiency2D", ";purity;efficiency", 110, 0., 1.1, 110, 0., 1.1);
  }
  if(fCheckClusters){
    fClusterPurity		= tfs->make<TH1D>("clusterPurity",     ";Purity;Clusters",     110, 0., 1.1);
    fClusterEfficiency		= tfs->make<TH1D>("clusterEfficiency", ";Efficiency;Clusters", 110, 0., 1.1);
    fClusterPurityEfficiency	= tfs->make<TH1D>("clusterPurityEfficiency", ";purityEfficiency;Clusters", 110, 0., 1.1);
    fClusterPurityEfficiency2D	= tfs->make<TH2D>("clusterPurityEfficiency2D", ";purity;efficiency", 110, 0., 1.1, 110, 0., 1.1);
  }

  fTree = tfs->make<TTree>("cheatertree","cheater tree");
  fTree->Branch("run",     &frun,     "run/I"    );
  fTree->Branch("event",   &fevent,   "event/I"  );
  fTree->Branch("trackid", &ftrackid, "trackid/I");
  fTree->Branch("pdg",     &fpdg,     "pdg/I"    );
  fTree->Branch("pmom",    &fpmom,    "pmom/D"   );
  fTree->Branch("hiteff",  &fhiteff,  "hiteff/D" );
  fTree->Branch("nclu",    &fnclu,    "nclu/I"   );
  fTree->Branch("clueff",  &fclueff              );
  fTree->Branch("clupur",  &fclupur              );
  fTree->Branch("cluid",   &fcluid               );
  fTree->Branch("nshw",    &fnshw,    "nshw/I"   );
  fTree->Branch("shweff",  &fshweff              );
  fTree->Branch("shwpur",  &fshwpur              );
  fTree->Branch("shwid",   &fshwid               );
  fTree->Branch("ntrk",    &fntrk,    "ntrk/I"   );
  fTree->Branch("trkeff",  &ftrkeff              );
  fTree->Branch("trkpur",  &ftrkpur              );
  fTree->Branch("trkid",   &ftrkid               );

  return;
}

//-------------------------------------------------------------------
// colID is the ID of the RecoBase object and colHits are the recob::Hits
// associated with it
void cheat::RecoCheckAna::CheckReco(int                                 const& colID,
				    std::vector< art::Ptr<recob::Hit> > const& allhits,
				    std::vector< art::Ptr<recob::Hit> > const& colHits,
				    std::map<std::pair<int, int>, std::pair<double, double> >& g4RecoBaseIDToPurityEfficiency)
{

  // grab the set of track IDs for these hits
  std::set<int> trackIDs = fBT->GetSetOfTrackIds(colHits);

  geo::View_t view = colHits[0]->View();

  std::set<int>::iterator itr = trackIDs.begin();
  while( itr != trackIDs.end() ){

    //std::cout << "*itr: " << *itr << std::endl;

    std::set<int> id;
    id.insert(*itr);

    // use the cheat::BackTrackerService to find purity and efficiency for these hits
    double purity     = fBT->HitCollectionPurity(id, colHits);
    double efficiency = fBT->HitCollectionEfficiency(id, colHits, allhits, view);

    // make the purity and efficiency pair
    std::pair<double, double> pe(purity, efficiency);

    // make the pair of the RecoBase object id to the pair of purity/efficiency
    std::pair<int, int> g4reco(*itr, colID);

    // insert idpe into the map
    g4RecoBaseIDToPurityEfficiency[g4reco] = pe;

    itr++;

  } // end loop over eveIDs

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::CheckRecoClusters(art::Event                                 const& evt,
					    std::string                                const& label,
					    art::Handle< std::vector<recob::Cluster> > const& clscol,
					    std::vector< art::Ptr<recob::Hit> >        const& allhits)
{

  art::FindManyP<recob::Hit> fmh(clscol, evt, label);

  for(size_t c = 0; c < clscol->size(); ++c){

    // get the hits associated with this event
    std::vector< art::Ptr< recob::Hit > > hits = fmh.at(c);

    this->CheckReco(clscol->at(c).ID(), allhits, hits, fG4ClusterIDToPurityEfficiency);

  }// end loop over clusters

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::CheckRecoTracks(art::Event                                 const& evt,
					  std::string                                const& label,
					  art::Handle< std::vector<recob::Track> >   const& tcol,
					  std::vector< art::Ptr<recob::Hit> >        const& allhits)
{

  art::FindManyP<recob::Hit> fmh(tcol, evt, label);

  for(size_t p = 0; p < tcol->size(); ++p){

    // get the hits associated with this event
    std::vector< art::Ptr< recob::Hit > > hits = fmh.at(p);

    this->CheckReco(tcol->at(p).ID(), allhits, hits, fG4TrackIDToPurityEfficiency);

  }// end loop over tracks

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::CheckRecoShowers(art::Event                                const& evt,
					   std::string                               const& label,
					   art::Handle< std::vector<recob::Shower> > const& scol,
					   std::vector< art::Ptr<recob::Hit> >       const& allhits)
{

  art::FindManyP<recob::Hit> fmh(scol, evt, label);

  for(size_t p = 0; p < scol->size(); ++p){

    // get the hits associated with this event
    std::vector< art::Ptr< recob::Hit > > hits = fmh.at(p);

    this->CheckReco(scol->at(p).ID(), allhits, hits, fG4ShowerIDToPurityEfficiency);

  }// end loop over events

  return;
}

//-------------------------------------------------------------------
//a true vertex will either consist of primary particles originating from
//the interaction vertex, or a primary particle decaying to make daughters
void cheat::RecoCheckAna::CheckRecoVertices(art::Event                                 const& evt,
					    std::string                                const& label,
					    art::Handle< std::vector<recob::Vertex> >  const& vtxcol,
					    std::vector< art::Ptr<recob::Hit> >        const& allhits)
{
  const sim::ParticleList& plist = fPI->ParticleList();

  std::vector< std::set<int> > ids(1);
  // loop over all primary particles and put their ids into the first set of the
  // vector.  add another set for each primary particle that also has daughters
  // and put those daughters into the new set
  // PartPair is a (track ID, particle pointer) pair
  for (const auto& PartPair: plist) {
    auto trackID = PartPair.first;
    if (!plist.IsPrimary(trackID)) continue;
    const simb::MCParticle& part = *(PartPair.second);
    ids[0].insert(trackID);
    if(part.NumberDaughters() > 0){
      std::set<int> dv;
      for(int d = 0; d < part.NumberDaughters(); ++d)
        dv.insert(part.Daughter(d));
      ids.push_back(std::move(dv));
    }//end if this primary particle has daughters
  }// end loop over primaries

  art::FindManyP<recob::Hit> fmh(vtxcol, evt, label);

  for(size_t v = 0; v < vtxcol->size(); ++v){

    // get the hits associated with this event
    std::vector< art::Ptr< recob::Hit > > hits = fmh.at(v);

    double maxPurity     = -1.;
    double maxEfficiency = -1.;

    for(size_t tv = 0; tv < ids.size(); ++tv){

      // use the cheat::BackTrackerService to find purity and efficiency for these hits
      double purity     = fBT->HitCollectionPurity(ids[tv], hits);
      double efficiency = fBT->HitCollectionEfficiency(ids[tv], hits, allhits, geo::k3D);

      if(purity     > maxPurity    ) maxPurity     = purity;
      if(efficiency > maxEfficiency) maxEfficiency = efficiency;
    }

    fVertexPurity    ->Fill(maxPurity);
    fVertexEfficiency->Fill(maxEfficiency);
    fVertexPurityEfficiency->Fill(maxPurity*maxEfficiency);

  }// end loop over vertices

  return;
}

//-------------------------------------------------------------------
// in this method one should loop over the primary particles from a given
// MCTruth collection
/// \todo need to divy it up in the case where there is more than 1 true interaction in a spill
void cheat::RecoCheckAna::CheckRecoEvents(art::Event                                 const& evt,
					  std::string                                const& label,
					  art::Handle< std::vector<recob::Event> >   const& evtcol,
					  std::vector< art::Ptr<recob::Hit> >        const& allhits)
{
  const sim::ParticleList& plist = fPI->ParticleList();

  // loop over all primaries in the plist and grab them and their daughters to put into
  // the set of track ids to pass on to the back tracker
  std::set<int> ids;
  for (const auto& PartPair: plist) {
    auto trackID = PartPair.first;
    if (!plist.IsPrimary(trackID)) continue;
    const simb::MCParticle& part = *(PartPair.second);
    ids.insert(trackID);
    for(int d = 0; d < part.NumberDaughters(); ++d)
      ids.insert(part.Daughter(d));
  }// end loop over primaries



  art::FindManyP<recob::Hit> fmh(evtcol, evt, label);

  for(size_t ev = 0; ev < evtcol->size(); ++ev){

    // get the hits associated with this event
    std::vector< art::Ptr< recob::Hit > > hits = fmh.at(ev);

    // use the cheat::BackTrackerService to find purity and efficiency for these hits
    double purity     = fBT->HitCollectionPurity(ids, hits);
    double efficiency = fBT->HitCollectionEfficiency(ids, hits, allhits, geo::k3D);

    fEventPurity    ->Fill(purity);
    fEventEfficiency->Fill(efficiency);
    fEventPurityEfficiency->Fill(purity*efficiency);

  }// end loop over events

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::FlattenMap(std::map<std::pair<int, int>, std::pair<double, double> > const& g4RecoBaseIDToPurityEfficiency,
				     std::map<int, std::vector<std::pair<int, std::pair<double, double> > > >& g4IDToRecoBasePurityEfficiency,
				     TH1D*                                                            purity,
				     TH1D*                                                            efficiency,
				     TH1D*                                                            purityEfficiency,
				     TH2D*                                                            purityEfficiency2D)
{

  std::map<std::pair<int, int>, std::pair<double, double> >::const_iterator rbItr = g4RecoBaseIDToPurityEfficiency.begin();

  // map of key cluster ID to pair of purity, efficiency
  std::map<int, std::pair<double, double> > recoBIDToPurityEfficiency;
  std::map<int, std::pair<double, double> >::iterator rbpeItr;

  while( rbItr != g4RecoBaseIDToPurityEfficiency.end() ){


    // trackID, cluster ID
    std::pair<int, int>       g4cl = rbItr->first;
    // purity, efficiency
    std::pair<double, double> pe   = rbItr->second;

    // add the efficiency and purity values for clusters corresponding
    // to the current g4 id to the map
    // pair of cluster id, pair of purity, efficiency
    std::pair<int, std::pair<double, double> > clpe(g4cl.second, pe);
    // g4IDToRecoBasePurityEfficiency is a map with key of trackID of a vector of clusterIDs of pairs of purity and efficiency
    g4IDToRecoBasePurityEfficiency[g4cl.first].push_back(clpe);

    // now find the maximum purity to determine the purity and efficiency
    // for this RecoBase object
    rbpeItr = recoBIDToPurityEfficiency.find(g4cl.second);
    if( rbpeItr != recoBIDToPurityEfficiency.end() ){
      std::pair<double, double> curpe = rbpeItr->second;
      if(pe.first > curpe.first) recoBIDToPurityEfficiency[g4cl.second] = pe;
    }
    else
      recoBIDToPurityEfficiency[g4cl.second] = pe;

    rbItr++;
  }

  rbpeItr = recoBIDToPurityEfficiency.begin();

  // now fill the histograms,
  while(rbpeItr != recoBIDToPurityEfficiency.end() ){
    purity    ->Fill(rbpeItr->second.first);
    efficiency->Fill(rbpeItr->second.second);
    purityEfficiency->Fill(rbpeItr->second.first*rbpeItr->second.second);
    purityEfficiency2D->Fill(rbpeItr->second.first,rbpeItr->second.second);
    rbpeItr++;
  }

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::FillResults(std::vector< art::Ptr<recob::Hit> > const& allhits)
{
  // map the g4 track id to energy deposited in a hit
  std::map<int, double> g4IDToHitEnergy;
  for(size_t h = 0; h < allhits.size(); ++h){
    const std::vector<sim::TrackIDE> hitTrackIDs = fBT->HitToTrackIDEs(allhits[h]);
    for(size_t e = 0; e < hitTrackIDs.size(); ++e){
      g4IDToHitEnergy[hitTrackIDs[e].trackID] += hitTrackIDs[e].energy;
    }
  } // end loop over hits to fill map

  // flatten the G4RecoBaseIDToPurityEfficiency maps to have just the g4ID as the key and the
  // rest of the information in vector form
  std::map<int, std::vector< std::pair<int, std::pair<double, double> > > > g4IDToClusterPurityEfficiency;
  std::map<int, std::vector< std::pair<int, std::pair<double, double> > > > g4IDToShowerPurityEfficiency;
  std::map<int, std::vector< std::pair<int, std::pair<double, double> > > > g4IDToTrackPurityEfficiency;
  std::map<int, std::vector< std::pair<int, std::pair<double, double> > > >::iterator g4peItr;

  if(fCheckClusters) this->FlattenMap(fG4ClusterIDToPurityEfficiency, g4IDToClusterPurityEfficiency, fClusterPurity, fClusterEfficiency, fClusterPurityEfficiency, fClusterPurityEfficiency2D);
  if(fCheckShowers)  this->FlattenMap(fG4ShowerIDToPurityEfficiency,  g4IDToShowerPurityEfficiency,  fShowerPurity,  fShowerEfficiency, fShowerPurityEfficiency, fShowerPurityEfficiency2D);
  if(fCheckTracks)   this->FlattenMap(fG4TrackIDToPurityEfficiency,   g4IDToTrackPurityEfficiency,   fTrackPurity,   fTrackEfficiency, fTrackPurityEfficiency, fTrackPurityEfficiency2D);

  // fill the tree vectors
  // get all the eveIDs from this event
  std::set<int> trackIDs = fBT->GetSetOfTrackIds();
  std::set<int>::const_iterator trackItr = trackIDs.begin();

  // loop over them
  while( trackItr != trackIDs.end() ){

    const simb::MCParticle* part = fPI->TrackIdToParticle_P(*trackItr);

    ftrackid = std::abs(*trackItr);
    fpdg     = part->PdgCode();
    fpmom    = part->P();

    // figure out how much of the energy deposited from this particle is stored in hits
    std::vector<const sim::IDE*> ides = fBT->TrackIdToSimIDEs_Ps(*trackItr);
    double totalDep = 0.;
    for(size_t i = 0; i < ides.size(); ++i) totalDep += ides[i]->energy;

    if(totalDep > 0.)
      fhiteff  = g4IDToHitEnergy[*trackItr]/totalDep;

    std::vector< std::pair<int, std::pair<double, double> > > clVec;
    std::vector< std::pair<int, std::pair<double, double> > > shVec;
    std::vector< std::pair<int, std::pair<double, double> > > trVec;


    if( g4IDToClusterPurityEfficiency.find(*trackItr) != g4IDToClusterPurityEfficiency.end() )
      clVec = g4IDToClusterPurityEfficiency.find(*trackItr)->second;

    if( g4IDToShowerPurityEfficiency.find(*trackItr) != g4IDToShowerPurityEfficiency.end() )
      shVec = g4IDToShowerPurityEfficiency.find(*trackItr)->second;

    if( g4IDToTrackPurityEfficiency.find(*trackItr) != g4IDToTrackPurityEfficiency.end() )
      trVec = g4IDToTrackPurityEfficiency.find(*trackItr)->second;

    fnclu = clVec.size();
    fnshw = shVec.size();
    fntrk = trVec.size();

    for(size_t c = 0; c < clVec.size(); ++c){
      fcluid .push_back(clVec[c].first);
      fclupur.push_back(clVec[c].second.first);
      fclueff.push_back(clVec[c].second.second);
    }

    for(size_t s = 0; s < shVec.size(); ++s){
      fshwid .push_back(shVec[s].first);
      fshwpur.push_back(shVec[s].second.first);
      fshweff.push_back(shVec[s].second.second);
    }

    for(size_t t = 0; t < trVec.size(); ++t){
      ftrkid .push_back(trVec[t].first);
      ftrkpur.push_back(trVec[t].second.first);
      ftrkeff.push_back(trVec[t].second.second);
    }

    fTree->Fill();

    trackItr++;
  }

  // clean up for the next event

  // clear the maps of G4 track id to efficiency and purity for
  // various RecoBase objects
  fG4ClusterIDToPurityEfficiency.clear();
  fG4ShowerIDToPurityEfficiency .clear();
  fG4TrackIDToPurityEfficiency  .clear();

  // clear the vectors hooked up to the tree
  fclueff.clear();
  fclupur.clear();
  fcluid .clear();
  ftrkeff.clear();
  ftrkpur.clear();
  ftrkid .clear();
  fshweff.clear();
  fshwpur.clear();
  fshwid .clear();

  return;
}


DEFINE_ART_MODULE(cheat::RecoCheckAna)
