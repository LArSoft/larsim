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


#include "larsim/LArG4/OpDetPhotonTable.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "larsim/Simulation/SimPhotons.h"

namespace larg4 {
  OpDetPhotonTable * TheOpDetPhotonTable;
  
  //--------------------------------------------------
  OpDetPhotonTable::OpDetPhotonTable()
  {
    fDetectedPhotons.clear();
  }

  //--------------------------------------------------
  OpDetPhotonTable * OpDetPhotonTable::Instance(bool /*LitePhotons*/ )
  {
    if(!TheOpDetPhotonTable){
      TheOpDetPhotonTable = new OpDetPhotonTable;
    }
    return TheOpDetPhotonTable;  
  }


  //--------------------------------------------------
  void OpDetPhotonTable::AddPhoton(size_t opchannel, sim::OnePhoton&& photon)
  {
    if( opchannel >= fDetectedPhotons.size() ) {

      std::cerr << "<<" << __PRETTY_FUNCTION__ << ">>"
		<< "\033[93m"
		<< "Invalid channel: " << opchannel
		<< "\033[00m"
		<< std::endl;
      throw std::exception();
    }
    fDetectedPhotons.at(opchannel).push_back(photon);
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
  void OpDetPhotonTable::ClearTable(const size_t nch)
  {
    if(fDetectedPhotons.size() != nch) fDetectedPhotons.resize(nch);
    for(size_t i=0; i<fDetectedPhotons.size(); ++i) {
      fDetectedPhotons.at(i).clear();
      fDetectedPhotons.at(i).SetChannel(i);
      //fDetectedPhotons.at(i).reserve(10000); // Just a guess on minimum # photons
    }

    for(std::map<int,std::map<int, int>>::iterator it=fLitePhotons.begin(); it!=fLitePhotons.end(); ++it)
      (it->second).clear();
    fLitePhotons.clear();
  }

  //--------------------------------------------------
  std::map<int, std::map<int, int>> OpDetPhotonTable::GetLitePhotons()
  {
    return fLitePhotons;
  }

  std::vector<sim::SimPhotons >& OpDetPhotonTable::GetPhotons() 
  { return fDetectedPhotons; }

  //--------------------------------------------------
  sim::SimPhotons& OpDetPhotonTable::GetPhotonsForOpChannel(size_t opchannel)
  {
    if(opchannel >= fDetectedPhotons.size()) {
      std::cerr << "<<" << __PRETTY_FUNCTION__ << ">>" 
		<< "Invalid channel Number: " << opchannel 
		<< std::endl;
    }
    return fDetectedPhotons.at(opchannel);
  }

  std::map<int,int>& OpDetPhotonTable::GetLitePhotonsForOpChannel(int opchannel)
  {
    return fLitePhotons[opchannel];
  }
  

}
