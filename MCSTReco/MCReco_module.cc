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
#include "art/Framework/Core/FindOneP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "MCShowerRecoAlg.h"
#include "MCTrackRecoAlg.h"

#include <memory>

class MCReco;

class MCReco : public art::EDProducer {
public:
  explicit MCReco(fhicl::ParameterSet const & p);
  virtual ~MCReco();

  void produce(art::Event & e) override;

private:

  // Declare member data here.
  std::string fG4Module;
  ::sim::MCRecoPart fPart;
  ::sim::MCRecoEdep fEdep;
  ::sim::MCShowerRecoAlg fMCSAlg;
  ::sim::MCTrackRecoAlg  fMCTAlg;
};

MCReco::MCReco(fhicl::ParameterSet const & pset)
  : fPart   (pset.get< fhicl::ParameterSet >("MCRecoPart"))
  , fEdep   (pset.get< fhicl::ParameterSet >("MCRecoEdep"))
  , fMCSAlg (pset.get< fhicl::ParameterSet >("MCShowerRecoAlg"))
  , fMCTAlg (pset.get< fhicl::ParameterSet >("MCTrackRecoAlg"))
{
  fG4Module = pset.get<std::string>("G4ModName","largeant");
  produces< std::vector< sim::MCShower> >();
  produces< std::vector< sim::MCTrack>  >();
  // Call appropriate produces<>() functions here.
}

MCReco::~MCReco()
{
  // Clean up dynamic memory and other resources here.
}

void MCReco::produce(art::Event & evt)
{
  std::unique_ptr< std::vector<sim::MCShower> > outShowerArray(new std::vector<sim::MCShower>);
  std::unique_ptr< std::vector<sim::MCTrack> > outTrackArray(new std::vector<sim::MCTrack>);

  // Retrieve mcparticles
  art::Handle<std::vector<simb::MCParticle> > mcpHandle;
  evt.getByLabel(fG4Module.c_str(),mcpHandle);
  if(!mcpHandle.isValid()) throw cet::exception(__FUNCTION__) << "Failed to retrieve simb::MCParticle";;

  // Find associations
  art::FindOneP<simb::MCTruth> ass(mcpHandle, evt, fG4Module.c_str());
  std::vector<simb::Origin_t> orig_array;
  orig_array.reserve(mcpHandle->size());
  for(size_t i=0; i<mcpHandle->size(); ++i) {
    const art::Ptr<simb::MCTruth> &mct = ass.at(i);
    orig_array.push_back(mct->Origin());
  }

  // Retrieve SimChannel
  art::Handle<std::vector<sim::SimChannel> > schHandle;
  evt.getByLabel(fG4Module.c_str(),schHandle);
  if(!schHandle.isValid()) throw cet::exception(__FUNCTION__) << "Failed to retrieve sim::SimChannel";

  const std::vector<simb::MCParticle>& mcp_array(*mcpHandle);
  fPart.AddParticles(mcp_array,orig_array);

  const std::vector<sim::SimChannel>&  sch_array(*schHandle);
  fEdep.MakeMCEdep(sch_array);

  fMCSAlg.Reconstruct(fPart,fEdep);

  fMCTAlg.Reconstruct(fPart,fEdep);

  for(auto const& mcs : fMCSAlg.MCShower())

    outShowerArray->push_back(mcs);

  for(auto const& mct : fMCTAlg.MCTrack())

    outTrackArray->push_back(mct);
    
  evt.put(std::move(outShowerArray));
  evt.put(std::move(outTrackArray));

}

DEFINE_ART_MODULE(MCReco)
