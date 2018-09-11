////////////////////////////////////////////////////////////////////////
// Class:       ShiftEdepSCE
// Plugin Type: producer (art v2_05_01)
// File:        ShiftEdepSCE_module.cc
//
// Generated at Thu Apr 19 00:41:18 2018 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "TNtuple.h"

namespace spacecharge {
  class ShiftEdepSCE;
}


class spacecharge::ShiftEdepSCE : public art::EDProducer {
public:
  explicit ShiftEdepSCE(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShiftEdepSCE(ShiftEdepSCE const &) = delete;
  ShiftEdepSCE(ShiftEdepSCE &&) = delete;
  ShiftEdepSCE & operator = (ShiftEdepSCE const &) = delete;
  ShiftEdepSCE & operator = (ShiftEdepSCE &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;
  void beginJob() override;

private:

  // Declare member data here.
  art::InputTag fEDepTag;
  bool          fMakeAnaTree;
  TNtuple*      fNtEdepAna;

    //IS calculationg
    larg4::ISCalcSeparate fISAlg;
};


spacecharge::ShiftEdepSCE::ShiftEdepSCE(fhicl::ParameterSet const & p)
  : fEDepTag(p.get<art::InputTag>("EDepTag")),
    fMakeAnaTree(p.get<bool>("MakeAnaTree",true))
{
  produces< std::vector<sim::SimEnergyDeposit> >();
}

void spacecharge::ShiftEdepSCE::beginJob()
{
  if(fMakeAnaTree){
    art::ServiceHandle<art::TFileService> tfs;
    fNtEdepAna = tfs->make<TNtuple>("nt_edep_ana","Edep PosDiff Ana Ntuple","energy:orig_x:orig_y:orig_z:orig_el:orig_ph:shift_x:shift_y:shift_z:shift_el:shift_ph");
  }

  art::ServiceHandle<sim::LArG4Parameters> lg4paramHandle;
  fISAlg.Initialize(lar::providerFrom<detinfo::LArPropertiesService>(),
		    lar::providerFrom<detinfo::DetectorPropertiesService>(),
		    &(*lg4paramHandle),
		    lar::providerFrom<spacecharge::SpaceChargeService>());

}

void spacecharge::ShiftEdepSCE::produce(art::Event & e)
{
  art::ServiceHandle<sim::LArG4Parameters> lg4paramHandle;
  fISAlg.Initialize(lar::providerFrom<detinfo::LArPropertiesService>(),
		    lar::providerFrom<detinfo::DetectorPropertiesService>(),
		    &(*lg4paramHandle),
		    lar::providerFrom<spacecharge::SpaceChargeService>());
  /*
  art::ServiceHandle<sim::LArG4Parameters> lg4paramHandle;
  fISAlg.Initialize(lar::providerFrom<detinfo::LArPropertiesService>(),
		    lar::providerFrom<detinfo::DetectorPropertiesService>(),
		    lg4paramHandle,
		    lar::providerFrom<spacecharge::SpaceChargeService>());
  */
  auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  std::unique_ptr< std::vector<sim::SimEnergyDeposit> > 
    outEdepVecPtr(new std::vector<sim::SimEnergyDeposit>() );
  auto & outEdepVec(*outEdepVecPtr);

  art::Handle< std::vector<sim::SimEnergyDeposit> > inEdepHandle;
  e.getByLabel(fEDepTag,inEdepHandle);
  auto const& inEdepVec(*inEdepHandle);
  
  outEdepVec.reserve(inEdepVec.size());
  
  geo::Vector_t posOffsetsStart{0.0,0.0,0.0};
  geo::Vector_t posOffsetsEnd{0.0,0.0,0.0};
  for(auto const& edep : inEdepVec){
    if(sce->EnableSimSpatialSCE()){
      posOffsetsStart = sce->GetPosOffsets({edep.StartX(),edep.StartY(),edep.StartZ()});
      posOffsetsEnd = sce->GetPosOffsets({edep.EndX(),edep.EndY(),edep.EndZ()});
    }
    fISAlg.Reset();
    fISAlg.CalculateIonizationAndScintillation(edep);
    outEdepVec.emplace_back(fISAlg.NumberScintillationPhotons(),
			    fISAlg.NumberIonizationElectrons(),
			    edep.Energy(),
			    geo::Point_t{(float)(edep.StartX()+posOffsetsStart.X()),
				(float)(edep.StartY()+posOffsetsStart.Y()),
				(float)(edep.StartZ()+posOffsetsStart.Z())},
			    geo::Point_t{(float)(edep.EndX()+posOffsetsEnd.X()),
				(float)(edep.EndY()+posOffsetsEnd.Y()),
				(float)(edep.EndZ()+posOffsetsEnd.Z())},
			    edep.StartT(),
			    edep.EndT(),
			    edep.TrackID(),
			    edep.PdgCode());
    if(fMakeAnaTree)
      fNtEdepAna->Fill(edep.Energy(),
		       edep.X(),edep.Y(),edep.Z(),edep.NumElectrons(),edep.NumPhotons(),
		       outEdepVec.back().X(),outEdepVec.back().Y(),outEdepVec.back().Z(),
		       outEdepVec.back().NumElectrons(),outEdepVec.back().NumPhotons());
  }

  e.put(std::move(outEdepVecPtr));
}

DEFINE_ART_MODULE(spacecharge::ShiftEdepSCE)
