////////////////////////////////////////////////////////////////////////
/// \file OpDetPhotonTable.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Implementation of the OpDetPhotonTable class.
//
// See comments in the OpDetPhotonTable.h file.
//
// Ben Jones, MIT, 11/12/12
//


#include "LArG4/OpDetPhotonTable.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/OpDetGeo.h"
#include "Simulation/SimPhotons.h"

namespace larg4 {
  OpDetPhotonTable * TheOpDetPhotonTable;
  
  //--------------------------------------------------
  OpDetPhotonTable::OpDetPhotonTable()
  {
    fDetectedPhotons.clear();
  }

  //--------------------------------------------------
  OpDetPhotonTable * OpDetPhotonTable::Instance(bool LitePhotons )
  {
    if(!TheOpDetPhotonTable){
      TheOpDetPhotonTable = new OpDetPhotonTable;
    }
    return TheOpDetPhotonTable;  
  }


  //--------------------------------------------------
  void OpDetPhotonTable::AddPhoton(int opchannel, sim::OnePhoton* photon)
  {
    if(!fDetectedPhotons[opchannel]) 
      {
	fDetectedPhotons[opchannel] = new sim::SimPhotons;
	fDetectedPhotons[opchannel]->SetChannel(opchannel);
      }
    fDetectedPhotons[opchannel]->push_back(*photon);
    //    mf::LogInfo("OpDetPhotonTable") << "Registering detection of a photon in opchannel " <<opchannel<<std::endl;
    
  }

  void OpDetPhotonTable::AddPhoton(std::map<int, std::map<int, int>>* StepPhotonTable)
  {
    for(auto it = StepPhotonTable->begin(); it!=StepPhotonTable->end(); it++)
    {
      for(auto in_it = it->second.begin(); in_it!=it->second.end(); in_it++)
      {
        fLitePhotons[it->first][in_it->first]+= in_it->second;
      }
    }
  }



  //--------------------------------------------------
  void OpDetPhotonTable::ClearTable()
  {
    for(std::map<int,sim::SimPhotons*>::iterator it=fDetectedPhotons.begin(); it!=fDetectedPhotons.end(); ++it)
      delete it->second;
    fDetectedPhotons.clear();
    for(std::map<int,std::map<int, int>>::iterator it=fLitePhotons.begin(); it!=fLitePhotons.end(); ++it)
      (it->second).clear();
    fLitePhotons.clear();
  }

  //--------------------------------------------------
  std::map<int, sim::SimPhotons* > OpDetPhotonTable::GetPhotons()
  {
    return fDetectedPhotons;
  }
  //--------------------------------------------------
  std::map<int, std::map<int, int>> OpDetPhotonTable::GetLitePhotons()
  {
    return fLitePhotons;
  }

  //--------------------------------------------------
  sim::SimPhotons* OpDetPhotonTable::GetPhotonsForOpChannel(int opchannel)
  {
    return fDetectedPhotons[opchannel];
  }
  std::map<int,int>& OpDetPhotonTable::GetLitePhotonsForOpChannel(int opchannel)
  {
    return fLitePhotons[opchannel];
  }
  

}
