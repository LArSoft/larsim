////////////////////////////////////////////////////////////////////////
/// \file VisualizationAction.h
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
/// VisualizationAction
/// 19-Mar-2002 Bill Seligman
///
/// Use UserAction to implement the standard visualization control for
/// a typical Geant4 job.  Everything in this class comes from the
/// Geant4 examples; the only difference is that it's put into an
/// UserAction class.
///
/// It uses multiple inheritance: it inherits from LArG4::UserAction,
/// in order to take advantage of Geant4's user hooks; it also
/// inherits from cfg::Observer, because it accesses a parameter from
/// an XML configuration file.
///
/// 04-Oct-2007 WGS: Accept a flag that controls the level of detail on
/// the event display.
///
/// 25-Feb-2009 WGS: Revised for FMWK/LArSoft

#ifndef LArG4_VisualizationAction_H
#define LArG4_VisualizationAction_H

#include "G4Base/UserAction.h"
#include "Geant4/globals.hh"

// Forward declarations of G4 classes that are arguments to
// UserAction methods.
class G4Run;
class G4Event;

namespace larg4 {

  class VisualizationAction : public g4b::UserAction
  {
  public:
    VisualizationAction();
    virtual ~VisualizationAction();

    virtual void BeginOfRunAction(const G4Run*); 
    virtual void EndOfRunAction(const G4Run*); 
    virtual void BeginOfEventAction(const G4Event*); 
    virtual void EndOfEventAction(const G4Event*); 

    /// Acessors, if we need them:
    G4double GetTrackEnergyCutoff() const {return m_energyCutoff;}
    void SetTrackEnergyCutoff(const G4double e) {m_energyCutoff = e;}

  private:
    /// Don't draw particles with energies less than this cut.
    G4double m_energyCutoff;

    /// Whether or not to draw neutral tracks (default is no).
    G4bool m_drawNeutrals;
  };

} // namespace LArG4

#endif // LArG4_VisualizationAction_H
