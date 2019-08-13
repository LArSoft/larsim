////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelID.h
/// \brief Unique identifier for a given LAr voxel.
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// This class defines a unique identifier for a given volume element
/// ("voxel") in the LAr volume.  It is sortable, can be tested for
/// equality, and is persistent under ROOT I/O.

/// (Actually, the term "voxel" is a mis-nomer, since we're also
/// keeping track of the time slice.  What's a four-dimensional volume
/// element?  A "tesseract element or "tessel"?)

#include "larsim/Simulation/LArVoxelID.h"
#include "larsim/Simulation/LArVoxelCalculator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <TLorentzVector.h>

#include <ostream>

namespace sim {

  //----------------------------------------------------------------------------
  /// Expert constructor based on actual bins.
  LArVoxelID::LArVoxelID( const int x,
			  const int y,
			  const int z,
			  const int t )
  {
    fbins = std::vector<int>(4);
    fbins[0] = x;
    fbins[1] = y;
    fbins[2] = z;
    fbins[3] = t;

  }

  //----------------------------------------------------------------------------
  /// Standard constructors.
  LArVoxelID::LArVoxelID( const TLorentzVector& coord )
  {
    // Copy each axis from the vector and convert it to a bin.
    fbins = std::vector<int>(4);
    art::ServiceHandle<sim::LArVoxelCalculator const> voxelCalc;
    for ( Ssiz_t a = 0; a != 4; ++a ){
      fbins[a] = voxelCalc->AxisToBin( a, coord[a] );
    }
  }

  //----------------------------------------------------------------------------
  LArVoxelID::LArVoxelID( const double x, const double y, const double z, const double t )
  {
    art::ServiceHandle<sim::LArVoxelCalculator const> voxelCalc;

    // Convert each axis into its corresponding bin.
    fbins = std::vector<int>(4);
    fbins[0] = voxelCalc->XAxisToBin( x );
    fbins[1] = voxelCalc->YAxisToBin( y );
    fbins[2] = voxelCalc->ZAxisToBin( z );
    fbins[3] = voxelCalc->TAxisToBin( t );
  }

  //----------------------------------------------------------------------------
  /// Destructor.
  LArVoxelID::~LArVoxelID() {}

  //----------------------------------------------------------------------------
  /// The accessors I expect to be used: The values of the
  /// co-ordinates at the bin centers.
  double LArVoxelID::X() const
  {
    art::ServiceHandle<sim::LArVoxelCalculator const> voxelCalc;
    return voxelCalc->XBinToAxis(fbins[0]);
  }

  //----------------------------------------------------------------------------
  double LArVoxelID::Y() const
  {
    art::ServiceHandle<sim::LArVoxelCalculator const> voxelCalc;
    return voxelCalc->YBinToAxis(fbins[1]);
  }

  //----------------------------------------------------------------------------
  double LArVoxelID::Z() const
  {
    art::ServiceHandle<sim::LArVoxelCalculator const> voxelCalc;
    return voxelCalc->ZBinToAxis(fbins[2]);
  }

  //----------------------------------------------------------------------------
  double LArVoxelID::T() const
  {
    art::ServiceHandle<sim::LArVoxelCalculator const> voxelCalc;
    return voxelCalc->TBinToAxis(fbins[3]);
  }

  //----------------------------------------------------------------------------
  double LArVoxelID::operator[]( const int i ) const
  {
    switch (i){
    case 0:
      return X(); break;
    case 1:
      return Y(); break;
    case 2:
      return Z(); break;
    case 3:
      return T(); break;
    }
    // I suppose I should put some error processing here; for now
    // I'll just return zero.
    return 0;
  }

  //----------------------------------------------------------------------------
  /// Put the contents on the output stream.  We have a choice: write
  /// the bin number, or write the position represented by the bins.
  /// For now, let's pick writing the positions.
  std::ostream& operator<< ( std::ostream& output, const LArVoxelID& id )
  {
    output << "(" << id.X()
	   << "," << id.Y()
	   << "," << id.Z()
	   << "," << id.T()
	   << ")";

    return output;
  }

  //----------------------------------------------------------------------------
  /// The comparison operator.  This a key function, since it
  /// establishes the sort order of the voxels in a list.
  bool LArVoxelID::operator<( const LArVoxelID& other ) const
  {
    // What is a good sort order for voxels in the list?  I'm not sure.
    // For now, pick an ordering but be prepared to change it: sort by
    // T, Z, X, then Y.

    if ( fbins[3] < other.fbins[3] ) return true;

    if ( fbins[3] == other.fbins[3] ){
      if ( fbins[2] < other.fbins[2] ) return true;

      if ( fbins[2] == other. fbins[2] ){
	if ( fbins[0] < other.fbins[0] ) return true;

	if ( fbins[0] == other.fbins[0] ){
	  if ( fbins[1] < other.fbins[1] ) return true;
	}
      }
    }

    return false;
  }

  //----------------------------------------------------------------------------
  /// Test for equality.  Handy, but not usually necessary.
  bool LArVoxelID::operator==( const LArVoxelID& other ) const
  {
    if ( fbins[0] == other.fbins[0]  &&
	 fbins[1] == other.fbins[1]  &&
	 fbins[2] == other.fbins[2]  &&
	 fbins[3] == other.fbins[3] )
      return true;

    return false;
  }

  //----------------------------------------------------------------------------
  LArVoxelID::operator TLorentzVector() const
  {
    return TLorentzVector( X(), Y(), Z(), T() );
  }

  //----------------------------------------------------------------------------
  LArVoxelID::operator TVector3() const
  {
    return TVector3( X(), Y(), Z() );
  }


} // namespace sim
