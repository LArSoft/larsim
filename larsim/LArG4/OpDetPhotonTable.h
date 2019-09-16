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
#ifndef OPDETPHOTONTABLE_h
#define OPDETPHOTONTABLE_h 1

#include <map>
#include <unordered_map>
#include <vector>
#include <string>

#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"

namespace larg4 {
  class OpDetPhotonTable
    {
    public:
      ~OpDetPhotonTable();
      static OpDetPhotonTable * Instance(bool LitePhotons = false);

      void AddPhoton( size_t opchannel, sim::OnePhoton&& photon, bool Reflected=false);
      void AddLitePhoton( int opchannel, int time, int nphotons, bool Reflected=false);
      void AddPhoton(std::map<int, std::map<int, int>>* StepPhotonTable, bool Reflected=false);
      void AddLitePhotons(std::map<int, std::map<int, int>>* StepPhotonTable, bool Reflected=false) { AddPhoton(StepPhotonTable, Reflected); }

      std::vector<sim::SimPhotons >& GetPhotons(bool Reflected=false) { return (Reflected ? fReflectedDetectedPhotons : fDetectedPhotons); }
      std::vector<sim::SimPhotons >& GetReflectedPhotons()            { return GetPhotons(true); }
      sim::SimPhotons&               GetPhotonsForOpChannel(size_t opchannel);
      sim::SimPhotons&               GetReflectedPhotonsForOpChannel(size_t opchannel);

      std::map<int, std::map<int, int> >    GetLitePhotons(bool Reflected=false) { return (Reflected ? fReflectedLitePhotons : fLitePhotons ); }
      std::map<int, std::map<int, int> >    GetReflectedLitePhotons()            { return GetLitePhotons(true); }
      std::map<int, int>&                   GetLitePhotonsForOpChannel(int opchannel)          { return fLitePhotons[opchannel]; }
      std::map<int, int>&                   GetReflectedLitePhotonsForOpChannel(int opchannel) { return fReflectedLitePhotons[opchannel]; }
      void ClearTable(size_t nch=0);

      void AddOpDetBacktrackerRecord(sim::OpDetBacktrackerRecord soc, bool Reflected=false);
    //  std::vector<sim::OpDetBacktrackerRecord>& GetOpDetBacktrackerRecords(); //Replaced by YieldOpDetBacktrackerRecords()
      std::vector<sim::OpDetBacktrackerRecord> YieldOpDetBacktrackerRecords();
      std::vector<sim::OpDetBacktrackerRecord> YieldReflectedOpDetBacktrackerRecords();


      void ClearEnergyDeposits();
      void AddEnergyDeposit(int n_elec,int n_photon,
			    double energy,
			    float start_x,float start_y, float start_z,
			    float end_x,float end_y,float end_z,
			    double start_time,double end_time,
			    int trackid,int pdgcode,
			    std::string const& vol="EMPTY");
      std::unordered_map<std::string, std::vector<sim::SimEnergyDeposit> > const& GetSimEnergyDeposits() const;
      //std::vector<sim::SimEnergyDeposit> & GetSimEnergyDeposits();

    protected:
      OpDetPhotonTable();

    private:

      void AddOpDetBacktrackerRecord(std::vector< sim::OpDetBacktrackerRecord > & RecordsCol,
                                     std::map<int, int> &ChannelMap,
                                     sim::OpDetBacktrackerRecord soc);


      std::map<int, std::map<int,int> >     fLitePhotons;
      std::map<int, std::map<int,int> >     fReflectedLitePhotons;
      std::vector< sim::OpDetBacktrackerRecord >      cOpDetBacktrackerRecordsCol; //analogous to scCol for electrons
      std::vector< sim::OpDetBacktrackerRecord >      cReflectedOpDetBacktrackerRecordsCol; //analogous to scCol for electrons
      std::map<int, int>  cOpChannelToSOCMap; //Where each OpChan is.
      std::map<int, int>  cReflectedOpChannelToSOCMap; //Where each OpChan is.
      std::vector<sim::SimPhotons> fDetectedPhotons;
      std::vector<sim::SimPhotons> fReflectedDetectedPhotons;


      std::unordered_map<std::string, std::vector<sim::SimEnergyDeposit> > fSimEDepCol;


    };

}


#endif
