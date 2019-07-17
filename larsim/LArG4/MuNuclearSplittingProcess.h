////////////////////////////////////////////////////////////////////////
/// \file  MuNuclearSplittingProcess.h
/// \brief Check of Geant4 to run the LArSoft detector simulation
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_MNSP_H
#define LARG4_MNSP_H

#include "Geant4/G4Types.hh"
#include "Geant4/G4WrapperProcess.hh"

class G4Step;
class G4Track;
class G4VParticleChange;

namespace larg4 {

class MuNuclearSplittingProcess : public G4WrapperProcess {
// Override PostStepDoIt method
  public:
    MuNuclearSplittingProcess() {};
    ~MuNuclearSplittingProcess() {};

    void SetNSplit(G4int nTrx) {fNSplit = nTrx;};
    void SetIsActive(G4bool doIt) {fActive = doIt;};

  private:
// Data members
    G4int fNSplit;
    G4bool fActive;
    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);

};


}// end namespace

#endif // MNSP
