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
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace larg4 {
  OpDetPhotonTable * TheOpDetPhotonTable;
  
  //--------------------------------------------------
  OpDetPhotonTable::OpDetPhotonTable()
  {
    fDetectedPhotons.clear();
  }
  OpDetPhotonTable::~OpDetPhotonTable(){}


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

  //--------------------------------------------------- cOpDetBacktrackerRecord population
  //J Stock. 11 Oct 2016
  void OpDetPhotonTable::AddOpDetBacktrackerRecord(sim::OpDetBacktrackerRecord soc){
    int iChan = soc.OpDetNum();
    //std::map<int,int> cOpChannelToSOCMap;
    std::map<int, int>::iterator channelPosition = cOpChannelToSOCMap.find(iChan);
    if (channelPosition == cOpChannelToSOCMap.end() ){
      cOpChannelToSOCMap[iChan] = cOpDetBacktrackerRecordsCol.size();
      cOpDetBacktrackerRecordsCol.emplace_back(std::move(soc));
    }else{
      unsigned int idtest = channelPosition->second;
      auto const& timePDclockSDPsMap = soc.timePDclockSDPsMap();
      for(auto const& timePDclockSDP : timePDclockSDPsMap){
        for(auto const& sdp : timePDclockSDP.second){
          double xyz[3] = {sdp.x, sdp.y, sdp.z};
          cOpDetBacktrackerRecordsCol.at(idtest).AddScintillationPhotons(
              sdp.trackID,
              timePDclockSDP.first,
              sdp.numPhotons,
              xyz,
              sdp.energy);
        }//end sdp : timesdp.second
      }//end const timesdp : timeSDPMap
    }// if chanPos == cOpChan else
    
    
    
  }//END void OpDetPhotonTable::AdOpDetBacktrackerRecords

  // cOpDetBacktrackerRecord return.
  std::vector<sim::OpDetBacktrackerRecord> OpDetPhotonTable::YieldOpDetBacktrackerRecords() {
    // we give the result to the caller, and don't retain it
    std::vector<sim::OpDetBacktrackerRecord> result;
    std::swap(result, cOpDetBacktrackerRecordsCol);
    cOpChannelToSOCMap.clear();
    return result;
  } // OpDetPhotonTable::YieldOpDetBacktrackerRecords()
  



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
  

  //--------------------------------------------------
  void OpDetPhotonTable::AddEnergyDeposit(int n_elec,int n_photon,
					  double energy,
					  float start_x,float start_y, float start_z,
					  float end_x,float end_y,float end_z,
					  double start_time,double end_time,
					  int trackid,int pdgcode,
					  std::string vol)
  {  
    fSimEDepCol[vol].emplace_back(n_elec,n_photon,
				  energy,
				  sim::SimEnergyDeposit::Point_t{start_x,start_y,start_z},
				  sim::SimEnergyDeposit::Point_t{end_x,end_y,end_z},
				  start_time,end_time,
				  trackid,pdgcode);
  }

  //--------------------------------------------------
  void OpDetPhotonTable::ClearEnergyDeposits()
  { fSimEDepCol.clear(); }
  
  
  //--------------------------------------------------
  std::unordered_map< std::string,std::vector<sim::SimEnergyDeposit> >& OpDetPhotonTable::GetSimEnergyDeposits()
  { return fSimEDepCol; }
  

}
