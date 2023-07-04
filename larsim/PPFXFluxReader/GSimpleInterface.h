#ifndef _GSIMPLEINTERFACE_H_
#define _GSIMPLEINTERFACE_H_

#include "larsim/PPFXFluxReader/FluxInterface.h"

#include "Tools/Flux/GNuMIFlux.h"
#include "Tools/Flux/GSimpleNtpFlux.h"

class TTree;
class TFile;

namespace fluxr {
  class GSimpleInterface : public FluxInterface {
  public:
    GSimpleInterface();
    ~GSimpleInterface();

    Long64_t GetEntries() const { return fNEntries; };
    int GetRun() const { return fRun; };
    float GetPOT() const { return fPOT; };
    TLorentzVector GetNuPosition() const { return fNuPos; };
    TLorentzVector GetNuMomentum() const { return fNuMom; };

    void SetRootFile(TFile* rootFileName);
    bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

  private:
    TTree* fFluxTree;
    TTree* fMetaTree;
    genie::flux::GSimpleNtpEntry* fGSimpleEntry;
    genie::flux::GSimpleNtpNuMI* fGSimpleNuMI;
    genie::flux::GSimpleNtpAux* fGSimpleAux;
    genie::flux::GSimpleNtpMeta* fGSimpleMeta;
    Long64_t fNEntries;
    int fRun;
    float fPOT;
    TLorentzVector fNuPos;
    TLorentzVector fNuMom;
  };

}

#endif // _GSIMPLEINTERFACE_H_
