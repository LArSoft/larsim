////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadout.h
/// \brief A Geant4 sensitive detector that accumulates voxel information.
///
/// \version $Id: LArVoxelReadout.h,v 1.2 2009/03/31 17:58:39 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// One way to implement voxels in Geant4 is to create a parallel
/// "read-out" geometry along with the real, physical geometry.  The
/// read-out geometry is implemented in LArVoxelReadoutGeometry; this
/// class is the sensitive detector for that geometry.  That is,
/// Geant4 will call this routine every time there is a step within a
/// volume of the read-out geometry; this routine then accumulates
/// information from that step.
///
/// In general, Geant4 expects to have per-event user information
/// attached to the G4Event in some way; their G4VSensitiveDetector
/// class supports this by allowing user-defined hit collections to
/// added to a G4HCOfThisEvent object (a collection of hit
/// collections; yes, it makes my head ache too!) that's part of each
/// G4Event.  
///
/// This class works differently, by accumulating the information in
/// its internal sim::LArVoxelList.  See LArVoxelListAction for how
/// this information is made available to the main LArG4 module.
/// 
/// Why define a parallel geometry?  Here are some reasons:
///
/// - The regular LAr TPC is one large volume of liquid argon.  When
///   Geant4 does its physics modeling, it can be unconstrained in
///   step size by the voxels.  Only for readout would the steps be
///   sub-divided.
///
/// - There may be more than one kind of readout, depending on a
///   detector's instrumentation (e.g., OpDets in addition to the wire
///   planes).  It's possible that the voxelization appropriate for
///   the wire planes may not be an appropriate readout for the other
///   readouts.  Geant4 allows the construction of multiple parallel
///   readouts, so this mechanism is relatively easy to extend for
///   each type of readout.

#ifndef LArG4_LArVoxelReadout_h
#define LArG4_LArVoxelReadout_h

#include <stdint.h>

#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/globals.hh"

#include "Simulation/SimChannel.h"
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/TimeService.h"
#include "Simulation/LArG4Parameters.h"
#include "LArG4/IonizationAndScintillation.h"



// Forward declarations
class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

namespace larg4 {

  class LArVoxelReadout : public G4VSensitiveDetector
  {
  public:
    // Constructor.  
    LArVoxelReadout(std::string const& name);

    // Destructor
    virtual ~LArVoxelReadout();

    // Required for classes that inherit from G4VSensitiveDetector.
    //
    // Called at start and end of each event.
    virtual void Initialize(G4HCofThisEvent*);
    virtual void EndOfEvent(G4HCofThisEvent*);

    // Called to clear any accumulated information.
    virtual void clear();

    // The key method of this class.  It's called by Geant4 for each
    // step within the read-out geometry.  It accumulates the energy
    // in the G4Step in the LArVoxelList.
    virtual G4bool ProcessHits( G4Step*, G4TouchableHistory* );

    // Empty methods; they have to be defined, but they're rarely
    // used in Geant4 applications.
    virtual void DrawAll();
    virtual void PrintAll();

    // Independent method; accumulates the information in an organized manner
    // call it after all MCTruth objects have been fed through Geant4
    void MakeSimChannelVector();

    // Independent method; clears the vector of SimChannels as well as the
    // channel number to SimChannel map.  Has to be separate from the 
    // clear method above because that run is run for every G4 event, ie
    // each MCTruth in the art::Event, while we want to only run this at 
    // the end of the G4 processing for each art::Event.
    void ClearSimChannelVector();

    // Independent method; returns the accumulated information
    const std::vector<sim::SimChannel>& GetSimChannels() const { return fChannels; }

  private:

    typedef enum radiologicaltype {
      notradiological,
      firstrad,
      subsequentrad } Radio_t;

    void DriftIonizationElectrons(G4ThreeVector stepMidPoint,
				  const double g4time,
				  int trackID,
				  Radio_t radiological=notradiological, 
				  unsigned int tickmax=4096); // used to randomize the TDC tick values  

    // Used in electron-cluster calculations
    // External parameters for the electron-cluster calculation.
    // obtained from LArG4Parameters, LArProperties, and DetectorProperties services
    double                                    fDriftVelocity[3];	       
    double 				      fLongitudinalDiffusion;  
    double 				      fTransverseDiffusion;    
    double 				      fElectronLifetime;       
    double 				      fElectronClusterSize;    
    double    				      fSampleRate; 	    
    double                                    fArgon39DecayRate;
    bool                                      fDontDriftThem;

    std::map<unsigned int, sim::SimChannel  > fChannelMap; ///< Map of channel number to SimChannel object
    std::vector<sim::SimChannel>              fChannels;   ///< Collection of SimChannels
    art::ServiceHandle<geo::Geometry>         fGeoHandle;  ///< Handle to the Geometry service
    art::ServiceHandle<sim::LArG4Parameters>  fLgpHandle;  ///< Handle to the LArG4 parameters service
    art::ServiceHandle<util::LArProperties>   fLarpHandle; ///< Handle to the LArProperties parameters service
    unsigned int                              fTPC;        ///< which TPC this LArVoxelReadout corresponds to
    unsigned int                              fCstat;      ///< and in which cryostat

    ::util::ElecClock                         fClock;      ///< TPC electronics clock
  };

}

#endif // LArG4_LArVoxelReadout_h
