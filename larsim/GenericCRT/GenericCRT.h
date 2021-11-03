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


#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"

#include <memory>

namespace sim{

    class GenericCRTUtility{

    public:

    unsigned int CopyAuxDetHittoAuxDetIDE(const sim::AuxDetHit &InputHit,sim::AuxDetIDE &outputIDE);

    unsigned int GetNumberofAuxDetChannels(const std::vector<sim::AuxDetHit> &InputHitVector,std::vector<unsigned int> &AuxDetChanNumbers );

    sim::AuxDetSimChannel const  GetAuxDetSimChannelByNumber(const std::vector<sim::AuxDetHit> &InputHitVector,unsigned int channel);

    void FillAuxDetSimChannels(const std::vector<sim::AuxDetHit> &InputHitVector, std::vector<sim::AuxDetSimChannel>  *AuxDetVector);

    private:

      art::ServiceHandle<geo::Geometry const> fGeo;

    };

}


#endif
