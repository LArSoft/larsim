////////////////////////////////////////////////////////////////////////
/// \file  G4BadIdeaAction.h
/// \brief this UserAction derived class is to implement catches to known bugs
///        in Geant4 that require grabbing const G4 objects and altering them -
///        a very bad idea in general.  Please do not add to this class without
///        discussing with the LArSoft Conveners
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

/// This class implements the LArG4::UserAction interface in order to
/// try to side step known bugs in version 4.9.4.p02 of Geant4
//
/// It uses multiple inheritance: it inherits from G4Base::UserAction,
/// in order to take advantage of Geant4's user hooks

#ifndef LArG4_G4BADIDEAACTION_H
#define LArG4_G4BADIDEAACTION_H

#include "nug4/G4Base/UserAction.h"

#include "Geant4/globals.hh"
#include <map>

#include "larcore/Geometry/Geometry.h"

// Forward declarations.
class G4Step;

namespace larg4 {

  class G4BadIdeaAction : public g4b::UserAction
  {
  public:
    // Standard constructors and destructors;
    G4BadIdeaAction(int );
    virtual ~G4BadIdeaAction();

    // UserActions method that we'll override, to obtain access to
    // Geant4's steps
    virtual void SteppingAction    (const G4Step*);

  private:

    art::ServiceHandle<geo::Geometry const> fGeo;  //< handle to geometry service
    int fNoIncomingMuons;

  };

} // namespace LArG4

#endif // LArG4_G4BADIDEAACTION_H
