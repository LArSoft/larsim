#include "larsim/LegacyLArG4/AllPhysicsLists.h"

namespace larg4 {
  AllPhysicsLists::AllPhysicsLists(detinfo::DetectorPropertiesData const& detProp)
    : fOpticalPhysics{detProp}
  {}
}
