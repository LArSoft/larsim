////////////////////////////////////////////////////////////////////////
// Class:       FilterSimPhotonTime
// Module Type: filter
// File:        FilterSimPhotonTime_module.cc
//
// Generated at Tue Jan 19 09:42:51 2016 by Wesley Ketchum using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include <memory>
#include <iostream>
#include <string>
#include <vector>

#include "lardataobj/Simulation/SimPhotons.h"

namespace simfilter {
  class FilterSimPhotonTime;
}

class simfilter::FilterSimPhotonTime : public art::SharedFilter {
public:
  explicit FilterSimPhotonTime(fhicl::ParameterSet const& p,
                               art::ProcessingFrame const&);

  // Plugins should not be copied or assigned.
  FilterSimPhotonTime(FilterSimPhotonTime const&) = delete;
  FilterSimPhotonTime(FilterSimPhotonTime&&) = delete;
  FilterSimPhotonTime& operator=(FilterSimPhotonTime const&) = delete;
  FilterSimPhotonTime& operator=(FilterSimPhotonTime&&) = delete;

private:
  bool filter(art::Event& e, art::ProcessingFrame const&) override;

  std::string const fSimPhotonsCollectionLabel;
  std::vector<std::vector<float>> const fTimeWindows;
  float const fMinTotalEnergy;
  float const fMinPhotonEnergy;
  bool const fDebug;
  std::size_t const fN;
  bool fUseReflectedPhotons;
  std::string fReflectedLabel;

  void CheckTimeWindows() const;
};

simfilter::FilterSimPhotonTime::FilterSimPhotonTime(
  fhicl::ParameterSet const& p,
  art::ProcessingFrame const&)
  : SharedFilter{p}
  , fSimPhotonsCollectionLabel(p.get<std::string>("SimPhotonsCollectionLabel"))
  , fTimeWindows(p.get<std::vector<std::vector<float>>>("TimeWindows"))
  , fMinTotalEnergy(p.get<float>("MinTotalEnergy", 0.0))
  , fMinPhotonEnergy(p.get<float>("MinPhotonEnergy", -1))
  , fDebug(p.get<bool>("Debug", false))
  , fN(fTimeWindows.size())
  , fUseReflectedPhotons(p.get<bool>("UseReflectedPhotons", false))
  , fReflectedLabel(p.get<std::string>("fReflectedLabel", "Reflected"))
{
  CheckTimeWindows();

  // For printing out debug messages, we want to serialize the
  // event-level calls so that the messages are not garbled.
  // Otherwise, this module works well for asynchronous event-level
  // calls.
  if (fDebug) {
    serialize();
  } else {
    async<art::InEvent>();
  }
}

void
simfilter::FilterSimPhotonTime::CheckTimeWindows() const
{

  if (fDebug)
    std::cout << "\tFilterSimPhotonTime: TimeWindows size is "
              << fTimeWindows.size() << std::endl;

  for (auto const& tw : fTimeWindows) {
    if (tw.size() != 2)
      throw cet::exception("FilterSimPhotonTime::CheckTimeWindows")
        << "Bad time window initialization: time window has wrong size (not 2)."
        << std::endl;

    if (fDebug)
      std::cout << "\t\tTimeWindow "
                << "[" << tw[0] << "," << tw[1] << "]" << std::endl;

    if (tw[0] > tw[1])
      throw cet::exception("FilterSimPhotonTime::CheckTimeWindows")
        << "Bad time window initialization: tw[0]>tw[1]. Reverse the order!"
        << std::endl;
  }
}

bool
simfilter::FilterSimPhotonTime::filter(art::Event& e,
                                       art::ProcessingFrame const&)
{
  auto const& simPhotonsCollection =
    *e.getValidHandle<std::vector<sim::SimPhotons>>(fSimPhotonsCollectionLabel);

  std::vector<double> sumEnergyArray(fN, 0.0);

  const std::vector<sim::SimPhotons> &simPhotonsCollectionReflected = fUseReflectedPhotons ?
    *e.getValidHandle<std::vector<sim::SimPhotons>>({fSimPhotonsCollectionLabel, fReflectedLabel}) : std::vector<sim::SimPhotons>();

  size_t n_sim_photons = simPhotonsCollection.size() + simPhotonsCollectionReflected.size();
    
  for (size_t i_pc = 0; i_pc < n_sim_photons; i_pc++) {
    const sim::SimPhotons &simphotons = (i_pc < simPhotonsCollection.size()) ? 
      simPhotonsCollection[i_pc] : simPhotonsCollectionReflected[i_pc - simPhotonsCollection.size()];

    if (fDebug)
      std::cout << "\tFilterSimPhotonTime: Processing simphotons channel "
                << simphotons.OpChannel() << std::endl;

    for (auto const& photon : simphotons)
      for (size_t i_tw = 0; i_tw < fN; ++i_tw) {
        auto const& tw(fTimeWindows[i_tw]);
        if (photon.Time >= tw[0] && photon.Time <= tw[1] &&
            photon.Energy > fMinPhotonEnergy) {

          if (fDebug) {
            std::string photon_string = (i_pc < simPhotonsCollection.size()) ? "Photon" : "Reflected Photon";
            std::cout << "\t\t" << photon_string << " with time " << photon.Time << " detected. "
                      << "Energy is  " << photon.Energy << "." << std::endl;
          }

          sumEnergyArray[i_tw] += photon.Energy;

          if (fDebug)
            std::cout << "\t\tTotal energy in this window (" << i_tw
                      << ") is now " << sumEnergyArray[i_tw] << std::endl;

          if (sumEnergyArray[i_tw] > fMinTotalEnergy)
            return true;
        }
      }
  }

  if (fDebug) {
    std::cout << "\tFilterSimPhotonTime: Final total energies are below min of "
              << fMinTotalEnergy << ":" << std::endl;
    for (size_t i_tw = 0; i_tw < fN; ++i_tw) {
      std::cout << "\t\tTimeWindow "
                << "[" << fTimeWindows[i_tw][0] << "," << fTimeWindows[i_tw][1]
                << "]: " << sumEnergyArray[i_tw] << std::endl;
    }
  }

  return false;
}

DEFINE_ART_MODULE(simfilter::FilterSimPhotonTime)
