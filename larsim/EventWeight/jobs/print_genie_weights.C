// Example gallery-based ROOT macro that retrieves and prints the contents of
// the MCEventWeight data product from an artroot-format file.
//
// Revised 17 March 2021 by Steven Gardiner <gardiner@fnal.gov>
//
// Make sure that the gallery ups product is set up before using this macro.
#include <iostream>
#include <string>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

const std::string genie_producer_label("generator");
const std::string eventweight_producer_label("genieeventweightTest");

void print_weights(const std::string& filename)
{

  std::vector<std::string> filenames{filename};

  gallery::Event ev(filenames);

  size_t event_count = 1;
  for (ev.toBegin(); !ev.atEnd(); ++ev) {

    std::cout << "art event " << event_count << '\n';

    // Use the GENIE producer label here
    auto mctruth_handle = ev.getValidHandle<std::vector<simb::MCTruth>>(genie_producer_label);
    auto gtruth_handle = ev.getValidHandle<std::vector<simb::GTruth>>(genie_producer_label);

    // Use the EventWeight producer label here
    auto weights_handle =
      ev.getValidHandle<std::vector<evwgh::MCEventWeight>>(eventweight_producer_label);

    // Loop through these objects for each neutrino vertex in the event
    for (size_t v = 0u; v < mctruth_handle->size(); ++v) {

      const simb::MCTruth& mc = mctruth_handle->at(v);
      const simb::GTruth& gt = gtruth_handle->at(v);
      const evwgh::MCEventWeight& mc_weights = weights_handle->at(v);

      // Extract desired information
      const auto& nu = mc.GetNeutrino();
      if (nu.Nu().NumberTrajectoryPoints() > 0) {
        double E_nu = nu.Nu().E(0);
        bool cc = (nu.CCNC() == 0);
        if (cc)
          std::cout << "CC";
        else
          std::cout << "NC";

        int mode = nu.Mode();
        if (mode == 0) {
          if (cc)
            std::cout << "QE";
          else
            std::cout << "EL";
        }
        else if (mode == 1)
          std::cout << "RES";
        else if (mode == 2)
          std::cout << "DIS";
        else if (mode == 3)
          std::cout << "COH";
        else if (mode == 10)
          std::cout << "MEC";
        else
          std::cout << " other";
        std::cout << " interaction with neutrino energy = " << E_nu << " GeV\n";
      }

      // Loop over all of the weights in the MCEventWeight object
      std::cout << "Weights\n";
      for (const auto& pair : mc_weights.fWeight) {
        std::string knob_name = pair.first;
        std::vector<double> weights = pair.second;

        std::cout << "  " << knob_name << ":\n";
        for (size_t u = 0u; u < weights.size(); ++u) {
          double w = weights.at(u);
          std::cout << "    universe #" << u << " has weight = " << w << '\n';
        }
      }
    }

    ++event_count;
  }
}
