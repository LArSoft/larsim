////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelCalculator.h
/// \brief Encapsulates calculation of LArVoxelID and LArVoxel parameters
/// \version $Id: ParticleHistory.cxx,v 1.1 2010/04/29 15:38:01 seligman Exp $///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// This class encapsulates the calculations associated with
/// computing the LArVoxelID, and provides access to the any LArVoxel
/// parameters from the input file(s).  

/// It is to be called using art::ServiceHandle<sim::LArVoxelCalculator> lvx;
/// The service makes it act like a singleton, but it knows about the
/// Parameters defined in the input file.

/// Definition: "Voxels" are three-dimensional "pixels"; basically
/// they divide the energy deposition in the LAr into (x,y,z) cubes.
/// Well, hyper-cubes actually, since we have to potentially include
/// divisions in time as well.

/// Example of a typical use:
///   const sim::LArVoxelCalculator* lhc = sim::LArVoxelCalculator::Instance();
///   double xSize = lhc->VoxelSizeX();

#ifndef sim_LArVoxelCalculator_H
#define sim_LArVoxelCalculator_H

#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace sim {

  class LArVoxelCalculator 
  {
  public:

    LArVoxelCalculator(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~LArVoxelCalculator();

    void reconfigure(fhicl::ParameterSet const& pset);

    /// Access to voxel dimensions and offsets.
    double VoxelSizeX()   const { return m_voxelSize[0];   }
    double VoxelSizeY()   const { return m_voxelSize[1];   }
    double VoxelSizeZ()   const { return m_voxelSize[2];   }
    double VoxelSizeT()   const { return m_voxelSize[3];   }
    double VoxelOffsetX() const { return m_voxelOffset[0]; }
    double VoxelOffsetY() const { return m_voxelOffset[1]; }
    double VoxelOffsetZ() const { return m_voxelOffset[2]; }
    double VoxelOffsetT() const { return m_voxelOffset[3]; }

    /// The energy in a voxel must be greater than this cut for it to
    /// be written to the output file.
    double EnergyCut() const { return m_energyCut; }

    /// Returns a step size that's reasonable to use so that each
    /// segment of a track will be contained within a single voxel.
    double SuggestedStepSize() const;

    /// Convert a co-ordinate axis (x, y, z, or t) into a bin number.
    /// The first argument is the axis (x=0, y=1, z=2, t=3) and the
    /// second is the value on that axis.
    int AxisToBin( const int, const double ) const;

    /// Provide an alternate access to the above routine with
    /// individual routines for the axes:
    int XAxisToBin( const double value ) const { return AxisToBin(0,value); }
    int YAxisToBin( const double value ) const { return AxisToBin(1,value); }
    int ZAxisToBin( const double value ) const { return AxisToBin(2,value); }
    int TAxisToBin( const double value ) const { return AxisToBin(3,value); }

    /// Get the value of an axis at the center of the given bin.  The
    /// first argument is the axis (x=0, y=1, z=2, t=3) and the second
    /// is the bin number on that axis.
    double BinToAxis( const int, const int ) const;

    /// Provide an alternate access to the above routine with
    /// individual routines for the axes:
    double XBinToAxis( const int value ) const { return BinToAxis(0,value); }
    double YBinToAxis( const int value ) const { return BinToAxis(1,value); }
    double ZBinToAxis( const int value ) const { return BinToAxis(2,value); }
    double TBinToAxis( const int value ) const { return BinToAxis(3,value); }

  private:

    typedef std::vector<double> vector_type;

    /// The sizes of the voxels in (x,y,z,t).  Units are (mm,ns).
    vector_type m_voxelSize;

    /// The offsets of the voxel binning from the origin in (x,y,z,t).
    /// Units are (mm,ns).
    vector_type m_voxelOffset;

    /// The total amount of energy in a voxel must be greater than
    /// this value for it to be written to the output.
    double m_energyCut;

  };

} // namespace sim

DECLARE_ART_SERVICE(sim::LArVoxelCalculator, LEGACY)
#endif //  sim_LArVoxelCalculator_H
