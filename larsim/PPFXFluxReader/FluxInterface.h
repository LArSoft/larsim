#ifndef _FLUXINTERFACE_H_
#define _FLUXINTERFACE_H_

#include "TLorentzVector.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace fluxr {
  class FluxInterface {
  public:
    virtual bool FillMCFlux(Long64_t ientry, simb::MCFlux& mclux) = 0;
    virtual const float GetPOT() = 0;
    virtual const Long64_t GetEntries() = 0;
    virtual const int GetRun() = 0;

    virtual const TLorentzVector GetNuPosition() = 0;
    virtual const TLorentzVector GetNuMomentum() = 0;
  };

}

#endif // _FLUXINTERFACE_H_
