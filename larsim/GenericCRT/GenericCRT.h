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

#include "larcorealg/Geometry/fwd.h"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/SimChannel.h"

#include <memory>
#include <string>
#include <vector>

namespace sim {

  class GenericCRTUtility {
  public:
    GenericCRTUtility(std::string const& energyUnitsScale,
                      geo::AuxDetGeometryCore const& auxDetGeom);

    AuxDetIDE toAuxDetIDE(AuxDetHit const& adhit) const;

    std::vector<unsigned int> GetAuxDetChannels(std::vector<AuxDetHit> const& adhits) const;

    std::vector<AuxDetSimChannel> GetAuxDetSimChannels(std::vector<AuxDetHit> const& adhits) const;

    AuxDetSimChannel GetAuxDetSimChannelByNumber(std::vector<AuxDetHit> const& adhits,
                                                 unsigned int inputchannel) const;

  private:
    geo::AuxDetGeometryCore const* fAuxDetGeom;
    double fEnergyUnitsScale;
  };

}

#endif
