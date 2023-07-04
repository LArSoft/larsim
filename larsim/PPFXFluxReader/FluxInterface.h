#ifndef _FLUXINTERFACE_H_
#define _FLUXINTERFACE_H_

#include "TLorentzVector.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace fluxr {
  class FluxInterface {
  public:
    virtual bool FillMCFlux(Long64_t ientry, simb::MCFlux& mclux) = 0;
    virtual float GetPOT() const = 0;
    virtual Long64_t GetEntries() const = 0;
    virtual int GetRun() const = 0;

    virtual TLorentzVector GetNuPosition() const = 0;
    virtual TLorentzVector GetNuMomentum() const = 0;
  };

}

#endif // _FLUXINTERFACE_H_
