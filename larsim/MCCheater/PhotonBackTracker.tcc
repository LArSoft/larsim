#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace cheat {

  //----------------------------------------------------------------
  template <typename Evt>
  bool PhotonBackTracker::CanRun(Evt const& evt)
  {
    return !evt.isRealData();
  }

  //----------------------------------------------------------------
  template <typename Evt>
  void PhotonBackTracker::PrepOpDetBTRs(Evt const& evt)
  {
    if (BTRsReady()) { return; }
    const std::vector<art::InputTag> G4ModuleLabels =
      (fG4ModuleLabels.empty()) ? std::vector{fG4ModuleLabel} : fG4ModuleLabels;

    auto compareBTRlambda = [](art::Ptr<sim::OpDetBacktrackerRecord> const& a,
                               art::Ptr<sim::OpDetBacktrackerRecord> const& b) {
      return a->OpDetNum() < b->OpDetNum();
    };

    for (auto& G4ModuleLabel : G4ModuleLabels) {
      auto const& btrHandle =
        evt.template getValidHandle<std::vector<sim::OpDetBacktrackerRecord>>(G4ModuleLabel);
      art::fill_ptr_vector(priv_OpDetBTRs, btrHandle);
      if (!std::is_sorted(priv_OpDetBTRs.begin(), priv_OpDetBTRs.end(), compareBTRlambda))
        std::sort(priv_OpDetBTRs.begin(), priv_OpDetBTRs.end(), compareBTRlambda);
    }
  }

  //----------------------------------------------------------------
  //ToDo: Figure out why I get OpHit* out of here instead of art::Ptr.
  template <typename Evt>
  void PhotonBackTracker::PrepOpFlashToOpHits(Evt const& evt)
  {
    if (OpFlashToOpHitsReady()) { return; }
    auto flashHandles = evt.template getMany<std::vector<recob::OpFlash>>();
    for (const auto& handle : flashHandles) {
      std::vector<art::Ptr<recob::OpFlash>> flash_vec;
      if (handle.failedToGet()) {
        mf::LogWarning("PhotonBackTracker")
          << " failed to get handle to recob::OpFlash. Has reco run yet?";
        return;
      }
      art::fill_ptr_vector(flash_vec, handle);
      auto tag = art::InputTag(handle.provenance()->moduleLabel());
      art::FindManyP<recob::OpHit> flash_hit_assn(flash_vec, evt, tag);
      for (size_t i = 0; i < flash_vec.size(); ++i) {
        art::Ptr<recob::OpFlash> flashp = flash_vec.at(i);
        std::vector<art::Ptr<recob::OpHit>> ophits = flash_hit_assn.at(i);
        auto check = priv_OpFlashToOpHits.emplace(flashp, ophits);
        if (check.second == false) {
          // loop ophit_ps
          // push_back to vector.
          for (auto& ophitp : ophits) {
            check.first->second.push_back(ophitp);
          }
        }
      }
    }
  }

  //----------------------------------------------------------------
  template <typename Evt>
  void PhotonBackTracker::PrepEvent(Evt const& evt)
  {
    if (!CanRun(evt)) {
      throw cet::exception("PhotonBackTracker") << "PhotonBackTracker cannot function."
                                                << "Is this file real data?";
    }
    priv_OpDetBTRs.clear();
    PrepOpDetBTRs(evt);
    PrepOpFlashToOpHits(evt);
  }
}
