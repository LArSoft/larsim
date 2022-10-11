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
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <memory>
#include <utility>

#include "TNtuple.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "larsim/Utils/SCEOffsetBounds.h"

namespace spacecharge {
  class ShiftEdepSCE;
}

class spacecharge::ShiftEdepSCE : public art::EDProducer {
public:
  explicit ShiftEdepSCE(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShiftEdepSCE(ShiftEdepSCE const&) = delete;
  ShiftEdepSCE(ShiftEdepSCE&&) = delete;
  ShiftEdepSCE& operator=(ShiftEdepSCE const&) = delete;
  ShiftEdepSCE& operator=(ShiftEdepSCE&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginJob() override;

private:
  // Declare member data here.
  art::InputTag fEDepTag;
  bool fMakeAnaTree;
  TNtuple* fNtEdepAna;

  //IS calculationg
  larg4::ISCalcSeparate fISAlg;
};

spacecharge::ShiftEdepSCE::ShiftEdepSCE(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fEDepTag(p.get<art::InputTag>("EDepTag"))
  , fMakeAnaTree(p.get<bool>("MakeAnaTree", true))
{
  produces<std::vector<sim::SimEnergyDeposit>>();
}

void spacecharge::ShiftEdepSCE::beginJob()
{
  if (fMakeAnaTree) {
    art::ServiceHandle<art::TFileService const> tfs;
    fNtEdepAna = tfs->make<TNtuple>(
      "nt_edep_ana",
      "Edep PosDiff Ana Ntuple",
      "energy:orig_x:orig_y:orig_z:orig_el:orig_ph:shift_x:shift_y:shift_z:shift_el:shift_ph");
  }
}

void spacecharge::ShiftEdepSCE::produce(art::Event& e)
{
  auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  auto const& inEdepVec = *e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEDepTag);

  auto outEdepVecPtr = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  auto& outEdepVec = *outEdepVecPtr;
  outEdepVec.reserve(inEdepVec.size());

  geo::Vector_t posOffsetsStart{0.0, 0.0, 0.0};
  geo::Vector_t posOffsetsEnd{0.0, 0.0, 0.0};
  for (auto const& edep : inEdepVec) {
    if (sce->EnableSimSpatialSCE()) {
      posOffsetsStart = sce->GetPosOffsets({edep.StartX(), edep.StartY(), edep.StartZ()});
      posOffsetsEnd = sce->GetPosOffsets({edep.EndX(), edep.EndY(), edep.EndZ()});
      if (larsim::Utils::SCE::out_of_bounds(posOffsetsStart) ||
          larsim::Utils::SCE::out_of_bounds(posOffsetsEnd)) {
        continue;
      }
    }
    auto const isData = fISAlg.CalcIonAndScint(detProp, edep);
    outEdepVec.emplace_back(
      isData.numPhotons,
      isData.numElectrons,
      0.0,
      edep.Energy(),
      geo::Point_t{(float)(edep.StartX() - posOffsetsStart.X()), //x should be subtracted
                   (float)(edep.StartY() + posOffsetsStart.Y()),
                   (float)(edep.StartZ() + posOffsetsStart.Z())},
      geo::Point_t{(float)(edep.EndX() - posOffsetsEnd.X()), //x should be subtracted
                   (float)(edep.EndY() + posOffsetsEnd.Y()),
                   (float)(edep.EndZ() + posOffsetsEnd.Z())},
      edep.StartT(),
      edep.EndT(),
      edep.TrackID(),
      edep.PdgCode());
    if (fMakeAnaTree)
      fNtEdepAna->Fill(edep.Energy(),
                       edep.X(),
                       edep.Y(),
                       edep.Z(),
                       edep.NumElectrons(),
                       edep.NumPhotons(),
                       outEdepVec.back().X(),
                       outEdepVec.back().Y(),
                       outEdepVec.back().Z(),
                       outEdepVec.back().NumElectrons(),
                       outEdepVec.back().NumPhotons());
  }

  e.put(std::move(outEdepVecPtr));
}

DEFINE_ART_MODULE(spacecharge::ShiftEdepSCE)
