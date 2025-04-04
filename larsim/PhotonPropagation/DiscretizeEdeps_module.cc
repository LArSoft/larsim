////////////////////////////////////////////////////////////////////////
// Class:       DiscretizeEdeps
// Plugin Type: producer (Unknown Unknown)
// File:        DiscretizeEdeps_module.cc
//
// Generated at Fri Apr  4 17:29:57 2025 by root using cetskelgen
// from cetlib version 3.18.02.
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
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include <memory>

namespace phot {
  class DiscretizeEdeps;
}


class phot::DiscretizeEdeps : public art::EDProducer {

using edep_tuple_t = std::tuple<int, int, int, int, int, int>;

public:
  explicit DiscretizeEdeps(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DiscretizeEdeps(DiscretizeEdeps const&) = delete;
  DiscretizeEdeps(DiscretizeEdeps&&) = delete;
  DiscretizeEdeps& operator=(DiscretizeEdeps const&) = delete;
  DiscretizeEdeps& operator=(DiscretizeEdeps&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fEDepLabel; // Label for the SimEnergyDeposit collection
  double fSpatialResolution; // Spatial resolution for discretization
  int fTimeResolution; // Time resolution for discretization Units: ns
};


phot::DiscretizeEdeps::DiscretizeEdeps(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fEDepLabel(p.get<art::InputTag>("EDepLabel", "largeant")),
    fSpatialResolution(p.get<double>("SpatialResolution", 1.0)), // Default: 1 cm
    fTimeResolution(p.get<double>("TimeResolution", 16)) // Default: 16 ns
{
  produces<std::vector<sim::SimEnergyDeposit>>();

}
struct EDepSum {
  int nsphotons = 0, nfphotons = 0, nelectrons = 0, pdg = 0;
  double energy = 0., syr = 1.;
};


void phot::DiscretizeEdeps::produce(art::Event& e) {

  //TODO -- CHECK ORIG TRACK ID
  auto output_edeps = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

  auto edeps = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEDepLabel);
  ///x, y, z, t, track ID
  std::map<edep_tuple_t, EDepSum> discretized_edeps;

  for (const auto& edep : *edeps) {
    // Discretize the energy deposit
    double x = edep.StartX();
    double y = edep.StartY();
    double z = edep.StartZ();
    double t = edep.StartT();

    // Discretize the position
    int discretized_x = std::floor(x / fSpatialResolution);
    int discretized_y = std::floor(y / fSpatialResolution);
    int discretized_z = std::floor(z / fSpatialResolution);

    // Discretize the time
    int discretized_t = std::floor(t / fTimeResolution);

    // Create a new SimEnergyDeposit with the discretized values
    edep_tuple_t the_tuple = std::make_tuple(discretized_x, discretized_y, discretized_z, discretized_t, edep.TrackID(), edep.OrigTrackID());
    if (discretized_edeps.find(the_tuple) == discretized_edeps.end()) {
      EDepSum the_edep_sum;
      the_edep_sum.pdg = edep.PdgCode();
      // the_edep_sum.x = discretized_x + 0.5 * fSpatialResolution;
      // the_edep_sum.y = discretized_y + 0.5 * fSpatialResolution;
      // the_edep_sum.z = discretized_z + 0.5 * fSpatialResolution;
      // the_edep_sum.t = discretized_t + 0.5 * fTimeResolution;
      discretized_edeps[the_tuple] = the_edep_sum;
    }
    EDepSum& sum = discretized_edeps[the_tuple];
    sum.nsphotons += edep.NumSPhotons();
    sum.nfphotons += edep.NumFPhotons();
    sum.nelectrons += edep.NumElectrons();
    sum.energy += edep.Energy();
  }

  for (const auto& [the_tuple, sum] : discretized_edeps) {
    int discretized_x = std::get<0>(the_tuple);
    int discretized_y = std::get<1>(the_tuple);
    int discretized_z = std::get<2>(the_tuple);
    int discretized_t = std::get<3>(the_tuple);
    int track_id = std::get<4>(the_tuple);
    int og_track_id = std::get<5>(the_tuple);
    double x = (discretized_x + .5) * fSpatialResolution;
    double y = (discretized_y + .5) * fSpatialResolution;
    double z = (discretized_z + .5) * fSpatialResolution;
    double t = (discretized_t + .5) * fTimeResolution;
    double edep = sum.energy;
    int nrelec = sum.nelectrons;
    int photons = (sum.nfphotons + sum.nsphotons);
    double syr = (photons > 0 ? 1.*sum.nfphotons/photons : 1.);
    // Add the new SimEnergyDeposit to the output collection

    output_edeps->push_back(
      sim::SimEnergyDeposit(
        photons, nrelec, syr/*1.*/, edep, {x, y, z}, {x, y, z}, t, t,
        track_id, sum.pdg, og_track_id
      )
    );
  }
  e.put(std::move(output_edeps));
}

DEFINE_ART_MODULE(phot::DiscretizeEdeps)
