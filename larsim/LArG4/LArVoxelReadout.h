////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadout.h
/// \brief A Geant4 sensitive detector that accumulates voxel information.
///
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
#include <vector>
#include <algorithm> // std::max()

#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/G4PVPlacement.hh"
#include "Geant4/globals.hh"

#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/LArG4/IonizationAndScintillation.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Forward declarations
class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;
namespace CLHEP { class HEPRandomEngine; }

namespace larg4 {

  /// Simple structure holding a TPC and cryostat number
  struct TPCID_t {
    unsigned short int Cryostat, TPC;
    bool operator< (const TPCID_t& than) const
      {
        return (Cryostat < than.Cryostat)
          || ((Cryostat == than.Cryostat) && (TPC < than.TPC));
      } // operator< ()
  }; // TPCID_t
  
  /**
   * @brief A G4PVPlacement with an additional identificator
   * @param IDTYPE type of ID class
   * 
   * This class is a G4PVPlacement with in addition an ID parameter.
   * The ID type is an object which can be default-constructed and copied,
   * better to be a POD.
   * 
   * This being a very stupid utility class, only the constructor that we
   * actually use is available. The others can be implemented in the same way.
   * Also the merry company of copy and move constuctors and operators is left
   * to the good will of the compiler, despite the destructor is specified.
   */
  template <class IDTYPE>
  class G4PVPlacementWithID: public G4PVPlacement {
      public:
    typedef IDTYPE ID_t;
    
    ID_t ID; ///< Physical Volume identificator
    
    /// Constructor
    G4PVPlacementWithID(const G4Transform3D& Transform3D, const G4String &pName,
      G4LogicalVolume* pLogical, G4VPhysicalVolume* pMother,
      G4bool pMany, G4int pCopyNo, G4bool pSurfChk = false,
      ID_t id = ID_t()
      ):
      G4PVPlacement
        (Transform3D, pName, pLogical, pMother, pMany, pCopyNo, pSurfChk),
      ID(id)
      {}
    
    /// Virtual destructor: does nothing more
    virtual ~G4PVPlacementWithID() {}
  }; // G4PVPlacementWithID<>
  
  /// A physical volume with a TPC ID
  typedef G4PVPlacementWithID<TPCID_t> G4PVPlacementInTPC;
  
  
  /**
   * @brief Transports energy depositions from GEANT4 to TPC channels.
   * 
   * This class acts on single energy depositions from GEANT4, simulating the
   * transportation of the ensuing ionisation electrons to the readout channels:
   * 
   * 1. the number of ionisation electrons is read from the current
   *   `larg4::IonizationAndScintillation` instance
   * 2. space charge displacement is optionally applied
   * 3. lifetime correction is applied
   * 4. charge is split in small electron clusters
   * 5. each cluster is subject to longitudinal and transverse diffusion
   * 6. each cluster is assigned to one TPC channel for each wire plane
   * 7. optionally, charge is forced to stay on the planes; otherwise charge
   *    drifting outside the plane is lost
   * 
   * For each energy deposition, entries on the appropriate `sim::SimChannel`
   * are added, with the information of the position where the energy deposit
   * happened (in global coordinates, centimeters), the ID of the Geant4
   * track which produced the deposition, and the quantized time of arrival to
   * the channel (in global TDC tick units).
   * At most one entry is added for each electron cluster, but entries from the
   * same energy deposit can be compacted if falling on the same TDC tick.
   * 
   * The main entry point of this class is the method `ProcessHits()`.
   * 
   * Options
   * --------
   * 
   * A few optional behaviours are supported:
   * 
   * * lead off-plane charge to the planes: regulated by
   *   `RecoverOffPlaneDeposit()`, if charge which reaches a wire plane
   *   is actually off it by less than the chosen margin, it's accounted for by
   *   that plane; by default the margin is 0 and all the charge off the plane
   *   is lost (with a warning)
   * 
   */
  class LArVoxelReadout : public G4VSensitiveDetector
  {
  public:
    /// Type of map channel -> sim::SimChannel
    typedef std::map<unsigned int, sim::SimChannel> ChannelMap_t;
    
    /// Collection of what it takes to set a `LArVoxelReadout` up.
    struct Setup_t {
      
      /// Random engine for charge propagation.
      CLHEP::HepRandomEngine* propGen = nullptr;
      
      /// Margin for charge recovery (see `LArVoxelReadout`).
      double offPlaneMargin = 0.0;
    }; // struct Setup_t
    
    
    /// Constructor. Can detect which TPC to cover by the name
    LArVoxelReadout(std::string const& name);
    
    /// Constructor. Sets which TPC to work on
    LArVoxelReadout
      (std::string const& name, unsigned int cryostat, unsigned int tpc);

    // Destructor
    virtual ~LArVoxelReadout();
    
    /// Reads all the configuration elements from `setupData`
    void Setup(Setup_t const& setupData);
    
    /// Associates this readout to one specific TPC
    void SetSingleTPC(unsigned int cryostat, unsigned int tpc);
    
    /// Sets this readout to discover the TPC of each processed hit
    void SetDiscoverTPC();
    
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

    // Independent method; clears the vector of SimChannels as well as the
    // channel number to SimChannel map.  Has to be separate from the 
    // clear method above because that run is run for every G4 event, ie
    // each MCTruth in the art::Event, while we want to only run this at 
    // the end of the G4 processing for each art::Event.
    void ClearSimChannels();

    /// Creates a list with the accumulated information for the single TPC
    std::vector<sim::SimChannel> GetSimChannels() const;
    
    /// Creates a list with the accumulated information for specified TPC
    std::vector<sim::SimChannel> GetSimChannels
      (unsigned short cryo, unsigned short tpc) const;

    //@{
    /// Returns the accumulated channel -> SimChannel map for the single TPC
    const ChannelMap_t& GetSimChannelMap() const;
    ChannelMap_t& GetSimChannelMap();
    //@}
    
    //@{
    /// Returns the accumulated channel -> SimChannel map for the specified TPC
    const ChannelMap_t& GetSimChannelMap
      (unsigned short cryo, unsigned short tpc) const;
    ChannelMap_t& GetSimChannelMap(unsigned short cryo, unsigned short tpc);
    //@}


  private:

    /**
     * @brief Sets the margin for recovery of charge drifted off-plane.
     * @param margin the extent of the margin on each frame coordinate [cm]
     * 
     * This method sets the margin for the recovery of off-plane ionization
     * charge. See `RecoverOffPlaneDeposit()` for a description of that feature.
     * 
     * This method is used by `LArVoxelReadout::Setup()`.
     */
    void SetOffPlaneChargeRecoveryMargin(double margin)
      { fOffPlaneMargin = std::max(margin, 0.0); }

    /// Sets the random generators to be used.
    void SetRandomEngines(CLHEP::HepRandomEngine* pPropGen);
    
    
    /**
     * @brief Returns the point on the specified plane closest to position.
     * @param pos the position to be tested (global coordinates, centimeters)
     * @param plane the plane to test the position against
     * @return a position on plane, unless pos is too far from it
     * 
     * This method considers the distance of the position `pos` from the active
     * part of the `plane` (see `geo::Plane::DeltaFromActivePlane()`).
     * If the position is less than a configurable margin far from the plane,
     * the closest point on the plane to that position is returned. Otherwise,
     * the position itself is returned.
     * 
     * Ionization charge may be drifted so that when it arrives to the plane, it
     * actually does not hit the area covered by wires. This can happen for many
     * reasons:
     * * space charge distortion led the point outside the fiducial volume
     *   (this may be prevented by specific code)
     * * diffusion pushes the charge outside the instrumented region
     * * the geometry of the wire planes is such that planes have different
     *   coverage and what one plane can cover, the next can't
     * 
     * The "recovery" consists in forcing the charge to the instrumented area
     * covered by the plane wires.
     * The distance of the drifted charge from each plane border is computed and
     * compared to the margin. If that distance is smaller than the margin, it
     * is neglected and the charge is assigned a new position on that border.
     * 
     * This method provides the position that should be used for the charge
     * deposition.
     * 
     * This is a simplistic approach to the simulation of border effects,
     * assuming that in fact the electric field, which is continuous and
     * pointing to the collection wires, will drive the charge to the wires even
     * when they are "off track".
     * No correction is applied for the additional time that such deviation
     * would take.
     * 
     */
    geo::Point_t RecoverOffPlaneDeposit
      (geo::Point_t const& pos, geo::PlaneGeo const& plane) const;
    
    void DriftIonizationElectrons(G4ThreeVector stepMidPoint,
                                  const double simTime,
                                  int trackID,
                                  unsigned short int cryostat, unsigned short int tpc);

    bool Has(std::vector<unsigned short int> v, unsigned short int tpc) const
    {  	
    	for (auto c: v) if (c == tpc) return true;	
    	return false;
    }

    // Used in electron-cluster calculations
    // External parameters for the electron-cluster calculation.
    // obtained from LArG4Parameters, LArProperties, and DetectorProperties services
    double                                    fDriftVelocity[3];
    double                                    fLongitudinalDiffusion;
    double                                    fTransverseDiffusion;
    double                                    fElectronLifetime;
    double                                    fElectronClusterSize;
    int					      fMinNumberOfElCluster;
    // for c2: unused private data members
    //double                                    fSampleRate;
    //int                                       fTriggerOffset;
    bool                                      fDontDriftThem;
    std::vector<unsigned short int>           fSkipWireSignalInTPCs;
    /// Charge deposited within this many [cm] from the plane is lead onto it.
    double                                    fOffPlaneMargin = 0.0;

    std::vector<std::vector<ChannelMap_t>>    fChannelMaps; ///< Maps of cryostat, tpc to channel data
    art::ServiceHandle<geo::Geometry>         fGeoHandle;  ///< Handle to the Geometry service
    art::ServiceHandle<sim::LArG4Parameters>  fLgpHandle;  ///< Handle to the LArG4 parameters service
    unsigned int                              fTPC;        ///< which TPC this LArVoxelReadout corresponds to
    unsigned int                              fCstat;      ///< and in which cryostat (if bSingleTPC is true)
    bool                                      bSingleTPC;  ///< true if this readout is associated with a single TPC
    
    CLHEP::HepRandomEngine*                   fPropGen = nullptr;  ///< random engine for charge propagation
    
    ::detinfo::ElecClock                         fClock;      ///< TPC electronics clock

    //these are the things for doing the separated EDeps
    void ProcessStep(G4Step*);
    
    G4ThreeVector                             fStepStart;
    G4ThreeVector                             fStepEnd;
    size_t fNSteps;
  };

}

#endif // LArG4_LArVoxelReadout_h
