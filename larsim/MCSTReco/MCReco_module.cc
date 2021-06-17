////////////////////////////////////////////////////////////////////////
// Class:       MCReco
// Module Type: producer
// File:        MCReco_module.cc
//
// Generated at Mon Aug 11 05:40:00 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "MCRecoEdep.h"
#include "MCShowerRecoAlg.h"
#include "MCTrackRecoAlg.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

class MCReco : public art::EDProducer {
public:
  explicit MCReco(fhicl::ParameterSet const & p);
//  virtual ~MCReco();

  void produce(art::Event & e) override;

private:

  // Declare member data here.
  art::InputTag fMCParticleLabel;
  art::InputTag fSimChannelLabel;
  bool fUseSimEnergyDeposit;

  ::sim::MCRecoPart fPart;
  ::sim::MCRecoEdep fEdep;
  ::sim::MCShowerRecoAlg fMCSAlg;
  ::sim::MCTrackRecoAlg  fMCTAlg;
};

MCReco::MCReco(fhicl::ParameterSet const & pset)
  : EDProducer{pset}
  , fPart   (pset.get< fhicl::ParameterSet >("MCRecoPart"))
  , fEdep   (pset.get< fhicl::ParameterSet >("MCRecoEdep"))
  , fMCSAlg (pset.get< fhicl::ParameterSet >("MCShowerRecoAlg"))
  , fMCTAlg (pset.get< fhicl::ParameterSet >("MCTrackRecoAlg"))
{

  //for backwards compatibility to using the "G4ModName" label...
  if(!(pset.get_if_present<art::InputTag>("MCParticleLabel",fMCParticleLabel) &&
       pset.get_if_present<art::InputTag>("SimChannelLabel",fSimChannelLabel)) ){

    mf::LogWarning("MCReco_module") << "USING DEPRECATED G4ModName CONFIG IN MCRECO_MODULE"
				    << "\nUse 'MCParticleLabel' and 'SimChannelLabel' instead.";

    fMCParticleLabel = pset.get<art::InputTag>("G4ModName","largeant");
    fSimChannelLabel = pset.get<art::InputTag>("G4ModName","largeant");
  }

  fUseSimEnergyDeposit = pset.get<bool>("UseSimEnergyDeposit",false);

  produces< std::vector< sim::MCShower> >();
  produces< std::vector< sim::MCTrack>  >();
  // Call appropriate produces<>() functions here.

//MCReco::~MCReco()
//{
  // Clean up dynamic memory and other resources here.
//}
}

void MCReco::produce(art::Event & evt)
{
//  std::unique_ptr< std::vector<sim::MCTrack> > outTrackArray(new std::vector<sim::MCTrack>);

  // Retrieve mcparticles
  art::Handle<std::vector<simb::MCParticle> > mcpHandle;
  evt.getByLabel(fMCParticleLabel,mcpHandle);
  if(!mcpHandle.isValid()) throw cet::exception(__FUNCTION__) << "Failed to retrieve simb::MCParticle";;

  // Find associations
  art::FindOneP<simb::MCTruth> ass(mcpHandle, evt, fMCParticleLabel);
  std::vector<simb::Origin_t> orig_array;
  orig_array.reserve(mcpHandle->size());
  for(size_t i=0; i<mcpHandle->size(); ++i) {
    const art::Ptr<simb::MCTruth> &mct = ass.at(i);
    orig_array.push_back(mct->Origin());
  }

  const std::vector<simb::MCParticle>& mcp_array(*mcpHandle);
  fPart.AddParticles(mcp_array,orig_array);


  // change implemented by David Caratelli to allow for MCRECO to run without SimChannels and using
  // SimEnergyDeposits instead
  if (fUseSimEnergyDeposit == true) {
    // Retrieve SimEnergyDeposit
    art::Handle<std::vector<sim::SimEnergyDeposit> > sedHandle;
    evt.getByLabel(fSimChannelLabel,sedHandle);
    if(!sedHandle.isValid()) throw cet::exception(__FUNCTION__) << "Failed to retrieve sim::SimEnergyDeposit";

    const std::vector<sim::SimEnergyDeposit>&  sed_array(*sedHandle);
    fEdep.MakeMCEdep(sed_array);
  }
  else {
    // Retrieve SimChannel
    art::Handle<std::vector<sim::SimChannel> > schHandle;
    evt.getByLabel(fSimChannelLabel,schHandle);
    if(!schHandle.isValid()) throw cet::exception(__FUNCTION__) << "Failed to retrieve sim::SimChannel";

    const std::vector<sim::SimChannel>&  sch_array(*schHandle);
    fEdep.MakeMCEdep(sch_array);
  }

  //Add MCShowers and MCTracks to the event
  evt.put(fMCSAlg.Reconstruct(fPart,fEdep));
  evt.put(fMCTAlg.Reconstruct(fPart,fEdep));

  fEdep.Clear();
  fPart.clear();
}

DEFINE_ART_MODULE(MCReco)
