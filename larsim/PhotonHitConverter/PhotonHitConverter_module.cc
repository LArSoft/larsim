////////////////////////////////////////////////////////////////////////
// Class:       PhotonHitConverter
// Plugin Type: producer (Unknown Unknown)
// File:        PhotonHitConverter_module.cc
//
// Generated at Tue Dec 13 05:30:41 2022 by Alejandro Castillo using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "artg4tk/pluginDetectors/gdml/PhotonHit.hh"

#include <memory>
#include <string>
#include <vector>

namespace sim {
  class PhotonHitConverter;
}

class sim::PhotonHitConverter : public art::EDProducer {
public:
  explicit PhotonHitConverter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  PhotonHitConverter(PhotonHitConverter const&) = delete;
  PhotonHitConverter(PhotonHitConverter&&) = delete;
  PhotonHitConverter& operator=(PhotonHitConverter const&) = delete;
  PhotonHitConverter& operator=(PhotonHitConverter&&) = delete;
  // Required functions.
  void produce(art::Event& e) override;

private:
  // Declare member data here.
  const bool fUseLitePhotons;
  size_t nOpChannels;
};

sim::PhotonHitConverter::PhotonHitConverter(fhicl::ParameterSet const& p)
  : EDProducer{p}, fUseLitePhotons(p.get<bool>("UseLitePhotons")) // ,
// More initializers here.
{
  // Call appropriate produces<>() functions here.
  if (fUseLitePhotons) {
    produces<std::vector<sim::SimPhotonsLite>>();
    produces<std::vector<sim::SimPhotonsLite>>("Reflected");
  }
  else {
    produces<std::vector<sim::SimPhotons>>();
    produces<std::vector<sim::SimPhotons>>("Reflected");
  }
}

void sim::PhotonHitConverter::produce(art::Event& e)
{
  geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
  nOpChannels = geom.NOpDets();

  //SimPhotonsLite
  std::unique_ptr<std::vector<sim::SimPhotonsLite>> photLiteCol{
    new std::vector<sim::SimPhotonsLite>{}};
  std::unique_ptr<std::vector<sim::SimPhotonsLite>> photLiteCol_ref{
    new std::vector<sim::SimPhotonsLite>{}};
  auto& photonLiteCollection{*photLiteCol};
  auto& photonLiteCollection_ref{*photLiteCol_ref};
  //SimPhotons
  std::unique_ptr<std::vector<sim::SimPhotons>> photCol{new std::vector<sim::SimPhotons>{}};
  std::unique_ptr<std::vector<sim::SimPhotons>> photCol_ref{new std::vector<sim::SimPhotons>{}};
  auto& photonCollection{*photCol};
  auto& photonCollection_ref{*photCol_ref};

  if (fUseLitePhotons) { //SimPhotonsLite
    photonLiteCollection.resize(nOpChannels);
    photonLiteCollection_ref.resize(nOpChannels);
    for (unsigned int i = 0; i < nOpChannels; ++i) {
      photonLiteCollection[i].OpChannel = i;
      photonLiteCollection_ref[i].OpChannel = i;
    }
  }
  else { //SimPhotons
    photonCollection.resize(nOpChannels);
    photonCollection_ref.resize(nOpChannels);
    for (unsigned int i = 0; i < nOpChannels; ++i) {
      photonCollection[i].fOpChannel = i;
      photonCollection_ref[i].fOpChannel = i;
    }
  }

  typedef std::vector<art::Handle<artg4tk::PhotonHitCollection>> HandleVector;
  auto allSims = e.getMany<artg4tk::PhotonHitCollection>();
  for (HandleVector::const_iterator i = allSims.begin(); i != allSims.end(); ++i) {
    const artg4tk::PhotonHitCollection& sims(**i);
    for (artg4tk::PhotonHitCollection::const_iterator j = sims.begin(); j != sims.end(); ++j) {
      const artg4tk::PhotonHit& photonHit = *j;
      if (fUseLitePhotons) {
        if (photonHit.GetEdep() > 6.19 * CLHEP::eV) {
          auto time = static_cast<int>(photonHit.GetTime());
          auto channel = static_cast<unsigned int>(photonHit.GetID());
          ++photonLiteCollection[channel].DetectedPhotons[time];
        }
        else {
          auto time = static_cast<int>(photonHit.GetTime());
          auto channel = static_cast<unsigned int>(photonHit.GetID());
          ++photonLiteCollection_ref[channel].DetectedPhotons[time];
        }
      }
      else {
        sim::OnePhoton photon;
        photon.SetInSD = false;
        photon.Energy = photonHit.GetEdep();
        auto time = photonHit.GetTime();
        photon.Time = time;
        auto channel = static_cast<unsigned int>(photonHit.GetID());
        if (photon.Energy > 6.19 * CLHEP::eV) {
          photonCollection[channel].insert(photonCollection[channel].end(), 1, photon);
        }
        else {
          photonCollection_ref[channel].insert(photonCollection_ref[channel].end(), 1, photon);
        }
      }
    }
  }
  if (fUseLitePhotons) {
    e.put(std::move(photLiteCol));
    e.put(std::move(photLiteCol_ref), "Reflected");
  }
  else {
    e.put(std::move(photCol));
    e.put(std::move(photCol_ref), "Reflected");
  }
}

DEFINE_ART_MODULE(sim::PhotonHitConverter)
