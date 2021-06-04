////////////////////////////////////////////////////////////////////////                                                           
//Class: ISTPC                                                                                                                     
//File: ISTPC.h and ISTPC.cxx                                                                                                      
//Description: Class including common functions needed for using the                                                               
//Hybrid Model for Photon Propagation                                                                                              
//UGR, 2021                                                                                                                        
//////////////////////////////////////////////////////////////////////// 



#include "larsim/IonizationScintillation/ISTPC.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace larg4 {
  //----------------------------------------------------------------------------                 
  ISTPC::ISTPC()
  {
    std::cout << "IonizationAndScintillation/ISTPC Initialize." << std::endl;
    std::cout << "Initializing the geometry of the detector." << std::endl;
    geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());

    fActiveVolumes = extractActiveVolumes(geom);
    {
      auto log = mf::LogTrace("IonAndScint") << "IonAndScint: active volume boundaries from "
                                             << fActiveVolumes.size() << " volumes:";
      for (auto const& [iCryo, box] : util::enumerate(fActiveVolumes)) {
        log << "\n - C:" << iCryo << ": " << box.Min() << " -- " << box.Max() << " cm";
      }
    } // local scope        
  }


  //---------------------------------------------------------------------------- 


 bool
 ISTPC::isScintInActiveVolume(geo::Point_t const& ScintPoint)
 {       
   return fActiveVolumes[0].ContainsPosition(ScintPoint);
 }
  //----------------------------------------------------------------------------                                 


                                                                                                                 
std::vector<geo::BoxBoundedGeo>
ISTPC::extractActiveVolumes(geo::GeometryCore const& geom)
{
  std::vector<geo::BoxBoundedGeo> activeVolumes;
  activeVolumes.reserve(geom.Ncryostats());

  for (geo::CryostatGeo const& cryo : geom.IterateCryostats()) {

    // can't use it default-constructed since it would always include origin                                                

    geo::BoxBoundedGeo box{cryo.TPC(0).ActiveBoundingBox()};

    for (geo::TPCGeo const& TPC : cryo.IterateTPCs())
      box.ExtendToInclude(TPC.ActiveBoundingBox());

    activeVolumes.push_back(std::move(box));

  } // for cryostats                                                                                                         

 return activeVolumes;
}

}
