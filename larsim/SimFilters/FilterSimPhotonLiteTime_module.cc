////////////////////////////////////////////////////////////////////////
// Class:       FilterSimPhotonLiteTime
// Module Type: filter
// File:        FilterSimPhotonLiteTime_module.cc
//
// Author: Gray Putnam -- ported from FilterSimPhotonTime
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
#include "lardataobj/Simulation/SimPhotons.h"

namespace simfilter {
  class FilterSimPhotonLiteTime;
}

class simfilter::FilterSimPhotonLiteTime : public art::SharedFilter {
public:
  explicit FilterSimPhotonLiteTime(fhicl::ParameterSet const& p,
                               art::ProcessingFrame const&);

  // Plugins should not be copied or assigned.
  FilterSimPhotonLiteTime(FilterSimPhotonLiteTime const&) = delete;
  FilterSimPhotonLiteTime(FilterSimPhotonLiteTime&&) = delete;
  FilterSimPhotonLiteTime& operator=(FilterSimPhotonLiteTime const&) = delete;
  FilterSimPhotonLiteTime& operator=(FilterSimPhotonLiteTime&&) = delete;

private:
  bool filter(art::Event& e, art::ProcessingFrame const&) override;

  std::string const fSimPhotonsLiteCollectionLabel;
  std::vector<std::vector<int>> const fTimeWindows;
  int const fMinTotalPhotons;
  bool const fDebug;
  std::size_t const fN;
  bool fUseReflectedPhotons;
  std::string fReflectedLabel;

  void CheckTimeWindows() const;
};

simfilter::FilterSimPhotonLiteTime::FilterSimPhotonLiteTime(
  fhicl::ParameterSet const& p,
  art::ProcessingFrame const&)
  : SharedFilter{p}
  , fSimPhotonsLiteCollectionLabel(p.get<std::string>("SimPhotonsLiteCollectionLabel"))
  , fTimeWindows(p.get<std::vector<std::vector<int>>>("TimeWindows"))
  , fMinTotalPhotons(p.get<int>("MinTotalPhotons"))
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
simfilter::FilterSimPhotonLiteTime::CheckTimeWindows() const
{

  if (fDebug)
    std::cout << "\tFilterSimPhotonLiteTime: TimeWindows size is "
              << fTimeWindows.size() << std::endl;

  for (auto const& tw : fTimeWindows) {
    if (tw.size() != 2)
      throw cet::exception("FilterSimPhotonLiteTime::CheckTimeWindows")
        << "Bad time window initialization: time window has wrong size (not 2)."
        << std::endl;

    if (fDebug)
      std::cout << "\t\tTimeWindow "
                << "[" << tw[0] << "," << tw[1] << "]" << std::endl;

    if (tw[0] > tw[1])
      throw cet::exception("FilterSimPhotonLiteTime::CheckTimeWindows")
        << "Bad time window initialization: tw[0]>tw[1]. Reverse the order!"
        << std::endl;
  }
}

bool
simfilter::FilterSimPhotonLiteTime::filter(art::Event& e,
                                       art::ProcessingFrame const&)
{
  auto const& simPhotonsLiteCollection =
    *e.getValidHandle<std::vector<sim::SimPhotonsLite>>(fSimPhotonsLiteCollectionLabel);

  std::vector<int> sumNPhotonArray(fN, 0);

  const std::vector<sim::SimPhotonsLite> &simPhotonsLiteCollectionReflected = fUseReflectedPhotons ?
    *e.getValidHandle<std::vector<sim::SimPhotonsLite>>({fSimPhotonsLiteCollectionLabel, fReflectedLabel}) : std::vector<sim::SimPhotonsLite>();

  size_t n_sim_photons = simPhotonsLiteCollection.size() + simPhotonsLiteCollectionReflected.size();

  if (fDebug) {
    std::cout << "New event to filter with total # sim photons: " << n_sim_photons << std::endl;
  }
    
  for (size_t i_pc = 0; i_pc < n_sim_photons; i_pc++) {
    const sim::SimPhotonsLite &simphotonslite = (i_pc < simPhotonsLiteCollection.size()) ? 
      simPhotonsLiteCollection[i_pc] : simPhotonsLiteCollectionReflected[i_pc - simPhotonsLiteCollection.size()];

    if (fDebug)
      std::cout << "\tFilterSimPhotonLiteTime: Processing simphotonslite channel "
                << simphotonslite.OpChannel << std::endl;

    for (auto const& photon_pair: simphotonslite.DetectedPhotons) {
      for (size_t i_tw = 0; i_tw < fN; ++i_tw) {
        auto const& tw(fTimeWindows[i_tw]);
        if (photon_pair.first >= tw[0] && photon_pair.first <= tw[1]) {

          if (fDebug) {
            std::string photon_string = (i_pc < simPhotonsLiteCollection.size()) ? "Photon" : "Reflected Photon";
            std::cout << "\t\t" << photon_string << " with number " << photon_pair.second 
                      << " at time " << photon_pair.first << " detected." << std::endl;
          }

          sumNPhotonArray[i_tw] += photon_pair.second;

          if (fDebug)
            std::cout << "\t\tTotal number of photons in this window (" << i_tw
                      << ") is now " << sumNPhotonArray[i_tw] << std::endl;

          if (sumNPhotonArray[i_tw] >= fMinTotalPhotons) 
            return true;
        }
      }
    }
  }

  if (fDebug) {
    std::cout << "\tFilterSimPhotonLiteTime: Final total numbers are below min of "
              << fMinTotalPhotons << ":" << std::endl;
    for (size_t i_tw = 0; i_tw < fN; ++i_tw) {
      std::cout << "\t\tTimeWindow "
                << "[" << fTimeWindows[i_tw][0] << "," << fTimeWindows[i_tw][1]
                << "]: " << sumNPhotonArray[i_tw] << std::endl;
    }
  }

  return false;
}

DEFINE_ART_MODULE(simfilter::FilterSimPhotonLiteTime)
