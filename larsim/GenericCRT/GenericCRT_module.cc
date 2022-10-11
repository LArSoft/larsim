////////////////////////////////////////////////////////////////////////
// Class:       GenericCRT
// Plugin Type: producer (art v3_05_01)
// File:        GenericCRT_module.cc
//
// Generated at Wed Oct 28 07:07:35 2020 by Andrzej Szelc using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "lardataobj/Simulation/AuxDetHit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "GenericCRT.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace sim {
  class GenericCRT;
}

class sim::GenericCRT : public art::EDProducer {
public:
  explicit GenericCRT(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GenericCRT(GenericCRT const&) = delete;
  GenericCRT(GenericCRT&&) = delete;
  GenericCRT& operator=(GenericCRT const&) = delete;
  GenericCRT& operator=(GenericCRT&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  std::string fEnergyUnitsScale;
  sim::GenericCRTUtility fCRTConvertUtil;

  // Declare member data here.
};

sim::GenericCRT::GenericCRT(fhicl::ParameterSet const& p)
  : EDProducer{p} //
  , fEnergyUnitsScale(p.get<std::string>("EnergyUnitsScale", "MeV"))
  , fCRTConvertUtil(fEnergyUnitsScale)
// More initializers here.
{

  produces<std::vector<sim::AuxDetSimChannel>>();
}

void sim::GenericCRT::produce(art::Event& e)
{
  // Implementation of required member function here.
  //std::unique_ptr< std::vector< sim::AuxDetSimChannel > > adCol (new  std::vector<sim::AuxDetSimChannel> );
  auto adCol = std::make_unique<std::vector<sim::AuxDetSimChannel>>();

  auto const& auxdethitcollection = e.getMany<std::vector<sim::AuxDetHit>>();

  for (size_t ii = 0; ii < auxdethitcollection.size(); ii++) {
    for (auto ch : fCRTConvertUtil.GetAuxDetSimChannels(*(auxdethitcollection.at(ii))))
      adCol->emplace_back(ch);
  }

  e.put(std::move(adCol));
}

DEFINE_ART_MODULE(sim::GenericCRT)
