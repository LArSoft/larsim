//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#ifndef LArStackingAction_H
#define LArStackingAction_H

#include "Geant4/globals.hh"
#include "Geant4/G4UserStackingAction.hh"
#include "Geant4/G4ThreeVector.hh"

class G4Track;

//#include "ExN04TrackerHit.hh"
//#include "ExN04MuonHit.hh"
class LArStackingActionMessenger;

class LArStackingAction : public G4UserStackingAction
{
  public:
    LArStackingAction(int );
    virtual ~LArStackingAction();

  public:
    // These 3 methods must be implemented by us. EC, 16-Feb-2011.
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();
    G4int fStack;

  private:
    //G4bool InsideRoI(const G4Track * aTrack,G4double ang);
    std::string InsideTPC(const G4Track * aTrack);
    //G4VHitsCollection* GetCollection(G4String colName);

    //ExN04TrackerHitsCollection* trkHits;
    //ExN04MuonHitsCollection* muonHits;
    //LArStackingActionMessenger* theMessenger;

    G4int fstage;
    G4int freqMuon;
    G4int freqIsoMuon;
    G4int freqIso;
    G4double fangRoI;



  public:
    inline void     SetNRequestMuon(G4int val)          { freqMuon = val;     }
    inline G4int    GetNRequestMuon()             const { return freqMuon;    }
    inline void     SetNRequestIsoMuon(G4int val)       { freqIsoMuon = val;  }
    inline G4int    GetNRequestIsoMuon()          const { return freqIsoMuon; }
    inline void     SetNIsolation(G4int val)            { freqIso = val;      }
    inline G4int    GetNIsolation()               const { return freqIso;     }
    inline void     SetRoIAngle(G4double val)           { fangRoI = val;      }
    inline G4double GetRoIAngle()                 const { return fangRoI;     }
};

#endif

