#ifndef _DK2NUINTERFACE_H_
#define _DK2NUINTERFACE_H_

#include "larsim/PPFXFluxReader/FluxInterface.h"
#include "dk2nu/tree/NuChoice.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

#include "fhiclcpp/ParameterSet.h"

class TTree;
class TFile;

#include "TLorentzRotation.h"
#include "TRandom3.h"

#include <vector>

namespace fluxr {
  class DK2NuInterface : public FluxInterface {
  public:
    Long64_t GetEntries() const { return fNEntries; };
    int GetRun() const { return fRun; };
    float GetPOT() const  { return fPOT; };
    TLorentzVector GetNuPosition() const { return fNuPos; };
    TLorentzVector GetNuMomentum() const { return fNuMom; };

    void SetRootFile(TFile* rootFile);
    bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

    bsim::Dk2Nu* GetDk2Nu() { return fDk2Nu; };
    bsim::NuChoice* GetNuChoice() { return fNuChoice; };

    void Init(fhicl::ParameterSet const& ps);
    void User2BeamPos(const TLorentzVector& usrxyz, TLorentzVector& beamxyz) const;
    void Beam2UserPos(const TLorentzVector& beamxyz, TLorentzVector& usrxyz) const;
    void Beam2UserP4(const TLorentzVector& beamp4, TLorentzVector& usrp4) const;
    TVector3 AnglesToAxis(double theta, double phi);

  private:
    TTree* fDk2NuTree;
    TTree* fDkMetaTree;
    bsim::Dk2Nu* fDk2Nu;
    bsim::DkMeta* fDkMeta;
    bsim::NuChoice* fNuChoice;
    Long64_t fNEntries;
    int fRun;
    float fPOT;

    TLorentzVector fNuPos;
    TLorentzVector fNuMom;

    TRotation fBeamRotXML, fTempRot;
    TLorentzRotation fBeamRot, fBeamRotInv;
    TVector3 fBeamPosXML;
    TLorentzVector fBeamZero;
    TVector3 detAV_rand_user;
    TLorentzVector fRandUser, fRandBeam;
  };

}

#endif // _DK2NUINTERFACE_H_
