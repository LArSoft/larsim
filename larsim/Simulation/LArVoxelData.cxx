////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelData.cxx
/// \brief Encapsulates the information we want store for a voxel.
///
/// \version $Id: LArVoxelData.cxx,v 1.2 2009/09/01 19:53:24 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "larsim/Simulation/LArVoxelData.h"

namespace sim {

  //----------------------------------------------------------------------------
  // Constructor; take care of any initializations.
  LArVoxelData::LArVoxelData() 
    : fenergy(0)
  {}
  
  //----------------------------------------------------------------------------
  // Destructor.
  LArVoxelData::~LArVoxelData() {}

  //----------------------------------------------------------------------------
  const LArVoxelData::key_type& LArVoxelData::TrackID( const size_type index ) const
  {
    const_iterator i = ftrackEnergy.begin();
    std::advance(i,index);
    return (*i).first;
  }

  //----------------------------------------------------------------------------
  const LArVoxelData::mapped_type& LArVoxelData::Energy( const size_type index ) const
  {
    const_iterator i = ftrackEnergy.begin();
    std::advance(i,index);
    return (*i).second;
  }

  //----------------------------------------------------------------------------
  void LArVoxelData::Add( const LArVoxelData& other )
  {
    // When we add one voxel's data to another, it becomes impossible
    // to keep the particle<->energy assignments anymore; the most
    // likely reason to add two voxels is because we're adding events
    // to make overlays, and so the particles' track numbers change.

    // So if we're adding another LArVoxelData to this one, move all
    // the energies to "unassigned" in the sum.

    this->RemoveAllTracks();
    fenergy += other.Energy();
  }

  //----------------------------------------------------------------------------
  LArVoxelData& LArVoxelData::operator*=( const double& value )
  {
    // Multiply all energies by the value.
    for ( iterator i = ftrackEnergy.begin(), end = ftrackEnergy.end(); i != end; ++i )
      {
	(*i).second *= value;
      }
    fenergy *= value;

    return (*this);
  }

  //----------------------------------------------------------------------------
  /// Just in case: define the result of "scalar * LArVoxelData" to be
  /// the same as "LArVoxelData * scalar".
  const LArVoxelData operator*(const double& value, const LArVoxelData& data) 
  {
    return LArVoxelData(data) *= value;
  }

  //----------------------------------------------------------------------------
  std::ostream& operator<< ( std::ostream& output, const LArVoxelData& data )
  {
    output << "Voxel: " << data.VoxelID() << std::endl;

    double unassigned = data.UnassignedEnergy();
    // Display the total energy then the breakdown of
    // the sum.
    output << data.Energy() << " = <ID,E>=";
    for ( LArVoxelData::const_iterator i = data.begin(); i != data.end(); ++i){
      if ( i != data.begin() )
	output << ",";
      
      output << "<" << (*i).first << "," << (*i).second << ">";
    }
    if( unassigned > 0 )
      output << ",<*," << unassigned << ">";

    return output;
  }

} // namespace sim
