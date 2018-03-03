////////////////////////////////////////////////////////////////////////
/// \file OpDetPhotonTable.h
//
/// \author  bjpjones@mit.edu
//  Eddited by JStock <jason.stock@mines.sdsmt.edu>
////////////////////////////////////////////////////////////////////////
// 
// This class holds a collection of PMT hits to be stored
// into an event as produced by the fast scintillation process.
//
// When a scintillating particle is stepped, it may generate
// one or more detected photons in each optical detector, 
// depending upon its location.
//
// For the fast scintillation process, the likelihood of generating
// a detected photo given a position in the detector is looked up.
// If one is detected, a hit in the relevant optical detector at
// an appropriate time (with ~20ns error bars for flight time) is 
// stored in this table, which is eventually read out at the end
// of the event by LArG4_module and stored in the event. 
//
// For slow scintillation / cerenkov processes, photons are
// generated and stepped about the detector, and if one steps
// into a volume designated at sensitive by Geant4, a simphoton
// is added to this table.
//
// The two sources can be distinguished by looking at the
// SetInSD flag of the OnePhoton object.
//
// Ben Jones, MIT, 11/10/12
//
//
//Changes have been made to this object to include the OpDetBacktrackerRecords for use in the photonbacktracker

#include "Geant4/G4PhysicalVolumeStore.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include <map>
#include <memory>
#include <exception>
#ifndef OPDETPHOTONTABLE_h
#define OPDETPHOTONTABLE_h 1

//#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace sim
{
  class OnePhoton; 
  class SimPhotons;
  class SimPhotonsLite;
  class SimEnergyDeposit;
}

namespace larg4 {
  class OpDetPhotonTable
    {
    public:
      ~OpDetPhotonTable();
      static OpDetPhotonTable * Instance(bool LitePhotons = false);
   
      void AddPhoton( size_t opchannel, sim::OnePhoton&& photon);
      void AddPhoton( std::map<int, std::map<int, int>>* StepPhoton);

      std::vector<sim::SimPhotons >& GetPhotons();
      sim::SimPhotons&               GetPhotonsForOpChannel(size_t opchannel);
      
      std::map<int, std::map<int, int> >   GetLitePhotons();
      std::map<int, int>&                  GetLitePhotonsForOpChannel(int opchannel); 
      void ClearTable(size_t nch=0);

      void AddOpDetBacktrackerRecord(sim::OpDetBacktrackerRecord soc);
    //  std::vector<sim::OpDetBacktrackerRecord>& GetOpDetBacktrackerRecords(); //Replaced by YieldOpDetBacktrackerRecords()
      std::vector<sim::OpDetBacktrackerRecord> YieldOpDetBacktrackerRecords();

      void ClearAndReserveEnergyDeposits(size_t reserve_size=0);
      void AddEnergyDeposit(int n_elec,int n_photon,
			    double energy,
			    float start_x,float start_y, float start_z,
			    float end_x,float end_y,float end_z,
			    double start_time,double end_time,
			    int trackid,int pdgcode);
      std::vector<sim::SimEnergyDeposit> & GetSimEnergyDeposits();
            
    protected:
      OpDetPhotonTable();

    private:

      std::map<int, std::map<int,int> >     fLitePhotons;
      std::vector< sim::OpDetBacktrackerRecord >      cOpDetBacktrackerRecordsCol; //analogous to scCol for electrons
      std::map<int, int>  cOpChannelToSOCMap; //Where each OpChan is.
      std::vector<sim::SimPhotons> fDetectedPhotons;


      std::vector<sim::SimEnergyDeposit>        fSimEDepCol;    


    };

}


#endif
