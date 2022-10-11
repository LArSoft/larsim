#ifndef GENERICCRTUTIL_H
#define GENERICCRTUTIL_H

/*!
 * Title:   GenericCRT Utility Class
 * Author:  Andrzej Szelc (andrzejs@fnal.gov)
 *
 * Description:
 * Class with Algorithms to convert sim::AuxDetHits to sim::AuxDetSimChannels
 *
 */

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/SimChannel.h"

#include <memory>
#include <string>
#include <vector>

namespace sim {

  class GenericCRTUtility {

  public:
    GenericCRTUtility(const std::string energyUnitsScale);

    sim::AuxDetIDE toAuxDetIDE(sim::AuxDetHit const& InputHit) const;

    std::vector<unsigned int> GetAuxDetChannels(
      const std::vector<sim::AuxDetHit>& InputHitVector) const;

    std::vector<sim::AuxDetSimChannel> GetAuxDetSimChannels(
      const std::vector<sim::AuxDetHit>& InputHitVector) const;

    sim::AuxDetSimChannel GetAuxDetSimChannelByNumber(
      const std::vector<sim::AuxDetHit>& InputHitVector,
      unsigned int inputchannel) const;

  private:
    art::ServiceHandle<geo::Geometry const> fGeo;

    double fEnergyUnitsScale;
  };

}

#endif
