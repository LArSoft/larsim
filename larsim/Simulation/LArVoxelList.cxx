////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelList.cxx
/// \brief Container of LArVoxelID, energy information.
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "larsim/Simulation/LArVoxelList.h"

#include <iterator>
#include <vector>
#include <iostream>
#include <cmath>

namespace sim {

  //----------------------------------------------------------------------------
  // Nothing special need be done for the constructor or destructor.
  LArVoxelList::LArVoxelList() {}

  //----------------------------------------------------------------------------
  LArVoxelList::~LArVoxelList() {}

  //----------------------------------------------------------------------------
  void LArVoxelList::Add( const LArVoxelList& other )
  {
    // Go through "other" list, adding its voxels to "this" list.
    for ( const_iterator i = other.m_voxelList.begin(); i != m_voxelList.end(); ++i ){
      m_voxelList[ (*i).first ] += (*i).second;
    }
  }

  //----------------------------------------------------------------------------
  LArVoxelList& LArVoxelList::operator*=( const double& value )
  {
    // Multiply each voxel energy by the value.
    for ( iterator i = m_voxelList.begin(); i != m_voxelList.end(); ++i ){
      (*i).second *= value;
    }
    return (*this);
  }

  //----------------------------------------------------------------------------
  /// Just in case: define the result of "scalar * LArVoxelList" to be
  /// the same as "LArVoxelList * scalar".
  const LArVoxelList operator*(const double& value, const LArVoxelList& list) 
  {
    return LArVoxelList(list) *= value;
  }

  //----------------------------------------------------------------------------
  /// Apply an energy cut to the voxels.
  void LArVoxelList::Cut( const double& cut )
  {
    // The safest way to do this is to create a list of voxel IDs that
    // fail the cut, then delete those IDs.

    // Define a list of IDs.
    typedef std::vector< key_type > keyList_type;
    keyList_type keyList;

    // Add each ID that fails the cut to the list.
    for ( const_iterator i = m_voxelList.begin(); i != m_voxelList.end(); ++i ){
      if ( (*i).second.Energy() < cut ) { 
	keyList.push_back( (*i).first ); 
      }
    }

    // Go through the list, deleting the voxels that are on the list.
    for ( keyList_type::const_iterator i = keyList.begin(); i != keyList.end(); ++i ){
      m_voxelList.erase( (*i) );
    }
  }

  //----------------------------------------------------------------------------
  const LArVoxelList::key_type& LArVoxelList::ID( const size_type index ) const
  {
    const_iterator i = m_voxelList.begin();
    std::advance(i,index);
    return (*i).first;
  }

  //----------------------------------------------------------------------------
  double LArVoxelList::Energy( const size_type index ) const
  {
    const_iterator i = m_voxelList.begin();
    std::advance(i,index);
    return (*i).second.Energy();
  }

  //----------------------------------------------------------------------------
  std::ostream& operator<< ( std::ostream& output, const LArVoxelList& list )
  {
    // Determine a field width for the voxel number.
    LArVoxelList::size_type numberOfVoxels = list.size();
    int numberOfDigits = (int) std::log10( (double) numberOfVoxels ) + 1;

    // A simple header.
    output.width( numberOfDigits );
    output << "#" << ": < ID, energy >" << std::endl; 

    // Write each voxel on a separate line.
    LArVoxelList::size_type nVoxel = 0;
    for ( LArVoxelList::const_iterator voxel = list.begin(); voxel != list.end(); ++voxel, ++nVoxel ){
      output.width( numberOfDigits );
      output << nVoxel << ": " 
	     << "< " << (*voxel).first 
	     << ", " << (*voxel).second 
	     << " >\n";
    }
    
    return output;
  }

} // namespace sim
