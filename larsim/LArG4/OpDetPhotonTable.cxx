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
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace larg4 {
  OpDetPhotonTable * TheOpDetPhotonTable;

  //--------------------------------------------------
  OpDetPhotonTable::OpDetPhotonTable()
  {
    fDetectedPhotons.clear();
    fReflectedDetectedPhotons.clear();
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
  void OpDetPhotonTable::AddPhoton(size_t opchannel, sim::OnePhoton&& photon, bool Reflected)
  {
    if( opchannel >= fDetectedPhotons.size() ) {

      std::cerr << "<<" << __PRETTY_FUNCTION__ << ">>"
		<< "\033[93m"
		<< "Invalid channel: " << opchannel
		<< "\033[00m"
		<< std::endl;
      throw std::exception();
    }
    if (!Reflected)
      fDetectedPhotons.at(opchannel).push_back(photon);
    else
      fReflectedDetectedPhotons.at(opchannel).push_back(photon);
  }

  //--------------------------------------------------
  void OpDetPhotonTable::AddLitePhoton( int opchannel, int time, int nphotons, bool Reflected)
  {
    if (!Reflected)
      fLitePhotons[opchannel][time] += nphotons;
    else
      fReflectedLitePhotons[opchannel][time] += nphotons;
  }

  //--------------------------------------------------
  void OpDetPhotonTable::AddPhoton(std::map<int, std::map<int, int>>* StepPhotonTable, bool Reflected)
  {
    for(auto it = StepPhotonTable->begin(); it!=StepPhotonTable->end(); it++)
    {
      for(auto in_it = it->second.begin(); in_it!=it->second.end(); in_it++)
      {
        if (!Reflected)
          fLitePhotons[it->first][in_it->first]+= in_it->second;
        else
          fReflectedLitePhotons[it->first][in_it->first]+= in_it->second;
      }
    }
  }

  //--------------------------------------------------- cOpDetBacktrackerRecord population
  //J Stock. 11 Oct 2016
  void OpDetPhotonTable::AddOpDetBacktrackerRecord(sim::OpDetBacktrackerRecord soc, bool Reflected){
//    std::cout << "DEBUG: Adding to " << (Reflected?"Reflected":"Direct") << " cOpDetBTR" << std::endl;
    if (!Reflected)
      AddOpDetBacktrackerRecord(cOpDetBacktrackerRecordsCol, cOpChannelToSOCMap, soc);
    else
      AddOpDetBacktrackerRecord(cReflectedOpDetBacktrackerRecordsCol, cReflectedOpChannelToSOCMap, soc);
  }

  //--------------------------------------------------- cOpDetBacktrackerRecord population
  void OpDetPhotonTable::AddOpDetBacktrackerRecord(std::vector< sim::OpDetBacktrackerRecord > & RecordsCol,
                                                   std::map<int, int> & ChannelMap,
                                                   sim::OpDetBacktrackerRecord soc) {
    int iChan = soc.OpDetNum();
    std::map<int, int>::iterator channelPosition = ChannelMap.find(iChan);
    if (channelPosition == ChannelMap.end() ){
      ChannelMap[iChan] = RecordsCol.size();
      RecordsCol.emplace_back(std::move(soc));
    }else{
      unsigned int idtest = channelPosition->second;
      auto const& timePDclockSDPsMap = soc.timePDclockSDPsMap();
      for(auto const& timePDclockSDP : timePDclockSDPsMap){
        for(auto const& sdp : timePDclockSDP.second){
          double xyz[3] = {sdp.x, sdp.y, sdp.z};
          RecordsCol.at(idtest).AddScintillationPhotons(
              sdp.trackID,
              timePDclockSDP.first,
              sdp.numPhotons,
              xyz,
              sdp.energy);
        }//end sdp : timesdp.second
      }//end const timesdp : timeSDPMap
    }// if chanPos == cOpChan else


//    std::cout << "DEBUG: Add to " << iChan << " to cOpDetBTR. Now " << RecordsCol.size() << " in size " << std::endl;
  }//END void OpDetPhotonTable::AdOpDetBacktrackerRecords


  //--------------------------------------------------
  // cOpDetBacktrackerRecord return.
  std::vector<sim::OpDetBacktrackerRecord> OpDetPhotonTable::YieldOpDetBacktrackerRecords() {
    // we give the result to the caller, and don't retain it
    std::vector<sim::OpDetBacktrackerRecord> result;
//    std::cout << "DEBUG: result.size()       = " << result.size() << std::endl;
//    std::cout << "DEBUG: cOpDetBTRCol.size() = " << cOpDetBacktrackerRecordsCol.size() << std::endl;
    std::swap(result, cOpDetBacktrackerRecordsCol);
//    std::cout << "DEBUG: std::swap(result, cOpDetBacktrackerRecordsCol);" << std::endl;
//    std::cout << "DEBUG: result.size()       = " << result.size() << std::endl;
//    std::cout << "DEBUG: cOpDetBTRCol.size() = " << cOpDetBacktrackerRecordsCol.size() << std::endl;
    cOpChannelToSOCMap.clear();
    return result;
  } // OpDetPhotonTable::YieldOpDetBacktrackerRecords()

  //--------------------------------------------------
  // cReflectedOpDetBacktrackerRecord return.
  std::vector<sim::OpDetBacktrackerRecord> OpDetPhotonTable::YieldReflectedOpDetBacktrackerRecords() {
    // we give the result to the caller, and don't retain it
    std::vector<sim::OpDetBacktrackerRecord> result;
//    std::cout << "DEBUG: result.size()           = " << result.size() << std::endl;
//    std::cout << "DEBUG: cReflOpDetBTRCol.size() = " << cReflectedOpDetBacktrackerRecordsCol.size() << std::endl;
    std::swap(result, cReflectedOpDetBacktrackerRecordsCol);
//    std::cout << "DEBUG: result.size()           = " << result.size() << std::endl;
//    std::cout << "DEBUG: cReflOpDetBTRCol.size() = " << cReflectedOpDetBacktrackerRecordsCol.size() << std::endl;
    cReflectedOpChannelToSOCMap.clear();
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
    if(fReflectedDetectedPhotons.size() != nch) fReflectedDetectedPhotons.resize(nch);
    for(size_t i=0; i<fReflectedDetectedPhotons.size(); ++i) {
      fReflectedDetectedPhotons.at(i).clear();
      fReflectedDetectedPhotons.at(i).SetChannel(i);
      //fDetectedPhotons.at(i).reserve(10000); // Just a guess on minimum # photons
    }

    for(auto it=fLitePhotons.begin(); it!=fLitePhotons.end(); ++it)
      (it->second).clear();
    for(auto it=fReflectedLitePhotons.begin(); it!=fReflectedLitePhotons.end(); ++it)
      (it->second).clear();
    fLitePhotons.clear();
    fReflectedLitePhotons.clear();
  }

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

  //--------------------------------------------------
  sim::SimPhotons& OpDetPhotonTable::GetReflectedPhotonsForOpChannel(size_t opchannel)
  {
    if(opchannel >= fReflectedDetectedPhotons.size()) {
      std::cerr << "<<" << __PRETTY_FUNCTION__ << ">>"
		<< "Invalid channel Number: " << opchannel
		<< std::endl;
    }
    return fReflectedDetectedPhotons.at(opchannel);
  }


  //--------------------------------------------------
  void OpDetPhotonTable::AddEnergyDeposit(int n_photon, int n_elec, double scint_yield, 
					  double energy,
					  float start_x,float start_y, float start_z,
					  float end_x,float end_y,float end_z,
					  double start_time,double end_time,
					  int trackid,int pdgcode,
					  std::string const& vol)
  {
    fSimEDepCol[vol].emplace_back(n_photon, n_elec, scint_yield,
				  energy,
				  geo::Point_t{start_x,start_y,start_z},
				  geo::Point_t{end_x,end_y,end_z},
				  start_time,end_time,
				  trackid,pdgcode);
  }

  //--------------------------------------------------
  void OpDetPhotonTable::ClearEnergyDeposits()
  { fSimEDepCol.clear(); }


  //--------------------------------------------------
  std::unordered_map< std::string,std::vector<sim::SimEnergyDeposit> > const& OpDetPhotonTable::GetSimEnergyDeposits() const
  { return fSimEDepCol; }


}
