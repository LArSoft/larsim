#ifndef larsim_LegacyLArG4_AllPhysicsLists_h
#define larsim_LegacyLArG4_AllPhysicsLists_h

// Register physics lists so that they can be loaded for use.  See the
// section "Using Custom Physics Modules" on the web page
// http://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/LArG4 for more
// details.

#include "larsim/LegacyLArG4/CustomPhysicsBuiltIns.hh"
#include "larsim/LegacyLArG4/CustomPhysicsLArSoft.h"
#include "larsim/LegacyLArG4/FastOpticalPhysics.h"
#include "larsim/LegacyLArG4/NeutronHPphysics.hh"
namespace detinfo {
  class DetectorPropertiesData;
}

namespace larg4 {
  class AllPhysicsLists {
  public:
    explicit AllPhysicsLists(detinfo::DetectorPropertiesData const& detProp);

  private:
    template <typename T>
    using Factory_t = CustomPhysicsFactory<T>;

    // G4 built-in lists
    Factory_t<G4ChargeExchangePhysics> fChargeExchange{"ChargeExchange"};
    Factory_t<G4DecayPhysics> fDecayPhysics{"Decay"};
    Factory_t<G4EmExtraPhysics> fSynchrotronAndGN{"SynchrotronAndGN"};
    Factory_t<G4EmLivermorePhysics> fLowEnergyEm{"LowEnergyEm"};
    Factory_t<G4EmPenelopePhysics> fpenelopeEmPhysics{"penelopeEm"};
    Factory_t<G4EmStandardPhysics> fEmPhysics{"Em"};
    Factory_t<G4HadronElasticPhysics> fHadronElasticPhysics{"HadronElastic"};
    Factory_t<G4HadronElasticPhysicsHP> fHadronElasticHPPhysics{"HadronElasticHP"};
    Factory_t<G4HadronElasticPhysicsPHP> fHadronElasticPHPPhysics{"HadronElasticPHP"};
    Factory_t<G4HadronPhysicsQGSP_BERT> fHadronPhysics{"Hadron"};
    Factory_t<G4HadronPhysicsQGSP_BERT_HP> fHadronHPPhysics{"HadronHP"};
    Factory_t<G4IonPhysics> fIonPhysics{"Ion"};
    Factory_t<G4NeutronTrackingCut> fNeutronTrackingCut{"NeutronTrackingCut"};
    Factory_t<G4StoppingPhysics> fStoppingPhysics{"Stopping"};

    // LArSoft lists
    Factory_t<FastOpticalPhysics> fFastOpticalPhysics{"FastOptical"};
    Factory_t<NeutronHPphysics> fNeutronHPPhysics{"NeutronHP"};
    Factory_t<OpticalPhysics> fOpticalPhysics;
  };
}

#endif
