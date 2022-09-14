////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadoutGeometry.h
/// \brief Define the "parallel" geometry that's seen by the LAr Voxels.
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// This class defines the parallel geometry that will be divided into
/// the three-dimensional voxels for the detector read-out.
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

#ifndef LArG4_LArVoxelReadoutGeometry_h
#define LArG4_LArVoxelReadoutGeometry_h

#include "Geant4/G4String.hh"
#include "Geant4/G4Transform3D.hh"
#include "Geant4/G4UserLimits.hh"
#include "Geant4/G4VUserParallelWorld.hh"
#include "larcore/Geometry/Geometry.h"
#include "larsim/LegacyLArG4/LArVoxelReadout.h"

// Forward declarations
class G4PhysicalVolume;

namespace larg4 {

  class LArVoxelReadoutGeometry : public G4VUserParallelWorld {
  public:
    /// Collection of all it takes to set up this object.
    struct Setup_t {

      /// Set up data for `LArVoxelReadout`.
      larg4::LArVoxelReadout::Setup_t readoutSetup;

    }; // struct Setup_t

    /// Constructor: sets up all its LArVoxelReadout instances.
    LArVoxelReadoutGeometry(const G4String& name, Setup_t const& setupData);

    /// N.B. Sets the relevant DetectorClocksData and
    /// DetectorPropertiesData for the event in flight.  See
    /// definition below.
    class Sentry;

    /// The key method in this class; creates a parallel world view of
    /// those volumes relevant to the LAr voxel readout.  Required of
    /// any class that inherits from G4VUserParallelWorld
    void Construct() override;

  private:
    G4VPhysicalVolume* FindNestedVolume(G4VPhysicalVolume* mother,
                                        G4Transform3D& motherTransform,
                                        G4Transform3D& daughterTransform,
                                        std::string const& daughterName,
                                        unsigned int expectedNum);

    art::ServiceHandle<geo::Geometry const> fGeo; ///< Handle to the geometry
    std::unique_ptr<G4UserLimits> fStepLimit;     ///< G4 doesn't handle memory management,
                                                  ///< so we have to

    /// Data for `LArVoxelReadout` setup.
    larg4::LArVoxelReadout* flarVoxelReadout{nullptr};
    larg4::LArVoxelReadout::Setup_t fReadoutSetupData;
  };

  // This is dirty, dirty.  Sigh.  Due to the limited interface of G4
  // and the prevalent use of services in LArSoft, there isn't a much
  // better alternative.
  class LArVoxelReadoutGeometry::Sentry {
  public:
    Sentry(LArVoxelReadoutGeometry* const readout,
           detinfo::DetectorClocksData const& clockData,
           detinfo::DetectorPropertiesData const& detProp) noexcept
      : readout_{readout}
    {
      readout_->flarVoxelReadout->SetClockData(&clockData);
      readout_->flarVoxelReadout->SetPropertiesData(&detProp);
    }
    ~Sentry() noexcept
    {
      readout_->flarVoxelReadout->SetClockData(nullptr);
      readout_->flarVoxelReadout->SetPropertiesData(nullptr);
    }

  private:
    LArVoxelReadoutGeometry* const readout_;
  };

} // namespace larg4

#endif // LArG4_LArVoxelReadoutGeometry_h
