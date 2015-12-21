////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelCalculator_plugin.cc
/// \brief Encapsulates calculation of LArVoxelID and LArVoxel parameters
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

/// This service encapsulates the calculations associated with
/// computing the LArVoxelID, and provides access to the any LArVoxel
/// parameters from the input file(s).

/// Definition: "Voxels" are three-dimensional "pixels"; basically
/// they divide the energy deposition in the LAr into (x,y,z) cubes.
/// Well, hyper-cubes actually, since we have to potentially include
/// divisions in time as well.

// Framework includes

#include "larsim/Simulation/LArVoxelCalculator.h"

#include <cmath>

namespace sim {

  //----------------------------------------------------------------------------
  LArVoxelCalculator::LArVoxelCalculator(fhicl::ParameterSet const& pset, 
					 art::ActivityRegistry & /* reg */) 
  {

    this->reconfigure(pset);

  }

  //----------------------------------------------------------------------------
  /// Destructor.
  LArVoxelCalculator::~LArVoxelCalculator()
  {
  }

  //----------------------------------------------------------------------------
  void LArVoxelCalculator::reconfigure(fhicl::ParameterSet const& pset)
  {
    m_energyCut = pset.get<double>("VoxelEnergyCut");

    double voxelSizeX   = pset.get<double>("VoxelSizeX");
    double voxelSizeY   = pset.get<double>("VoxelSizeY");
    double voxelSizeZ   = pset.get<double>("VoxelSizeZ");
    double voxelSizeT   = pset.get<double>("VoxelSizeT");
    double voxelOffsetX = pset.get<double>("VoxelOffsetX");
    double voxelOffsetY = pset.get<double>("VoxelOffsetY");
    double voxelOffsetZ = pset.get<double>("VoxelOffsetZ");
    double voxelOffsetT = pset.get<double>("VoxelOffsetT");
    
    m_voxelSize.clear();
    m_voxelSize.push_back( voxelSizeX );
    m_voxelSize.push_back( voxelSizeY );
    m_voxelSize.push_back( voxelSizeZ );
    m_voxelSize.push_back( voxelSizeT );
    m_voxelOffset.clear();
    m_voxelOffset.push_back( voxelOffsetX );
    m_voxelOffset.push_back( voxelOffsetY );
    m_voxelOffset.push_back( voxelOffsetZ );
    m_voxelOffset.push_back( voxelOffsetT );

    return;
  }

  //----------------------------------------------------------------------------
  /// Returns a Monte-Carlo step size that's reasonable to use so that
  /// each segment of a track will be contained within a single voxel.
  double LArVoxelCalculator::SuggestedStepSize() const
  {
    return std::min( m_voxelSize[0], std::min( m_voxelSize[1], m_voxelSize[2] ) );
  }

  //----------------------------------------------------------------------------
  /// Convert a co-ordinate axis (x, y, z, or t) into a bin number.
  /// The first argument is the axis (x=0, y=1, z=2, t=3) and the
  /// second is the value on that axis.
  int LArVoxelCalculator::AxisToBin( const int axis, 
				     const double coord ) const
  {
    // We have to be careful of how to handle the case when coord -
    // offset < 0.  The standard floor() function rounds the number in
    // the correct direction.
    return static_cast<int>( floor( ( coord - m_voxelOffset[axis] ) / m_voxelSize[axis] ) );
  }

  //----------------------------------------------------------------------------
  /// Get the value of an axis at the center of the given bin.  The
  /// first argument is the axis (x=0, y=1, z=2, t=3) and the second
  /// is the bin number on that axis.
  double LArVoxelCalculator::BinToAxis( const int axis, 
					const int bin ) const
  {
    return ( static_cast<double>(bin) + 0.5 ) * m_voxelSize[axis] + m_voxelOffset[axis];
  }

} // namespace sim
 
namespace sim {

  DEFINE_ART_SERVICE(LArVoxelCalculator)


} // namespace sim
 
