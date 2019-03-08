////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelID.h
/// \brief Unique identifier for a given LAr voxel.
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// This class defines a unique identifier for a given volume element
/// ("voxel") in the LAr volume.  It is sortable, can be tested for
/// equality, and is persistent under ROOT I/O.
///
/// (Actually, the term "voxel" is a mis-nomer, since we're also
/// keeping track of the time slice.  What's a four-dimensional volume
/// element?  A "tesseract element" or "tessel"?)
///
/// A LArVoxelID is created by supplying the (x,y,z,t) co-ordinates,
/// which can be in the form of a TLorentzVector:
///
///    sim::LArVoxelID id(1,2,3,4);
///    TLorentzVector v(1,2,3,4);
///    sim::LArVoxelID otherID(v);
///
/// There are several ways to "get" at the contents of a LArVoxelID:
///
///    sim::LArVoxelID id;
///    double x = id.X();           // axis-by-axis
///    int xBinNumber = id.XBin();  // bin number; might be useful for debugging
///    TLorentzVector pos(id);      // convert to TLorentzVector
///
/// Note that the first and third methods above both return the
/// bin-center value(s).
///
/// If you need to know the sizes of the bins (voxels), see
/// LArVoxelCalculator.
///
/// Output methods are implemented so this class can be a good
/// "citizen":
///
///    sim::LArVoxelID voxelID;
///    std::cout << voxelID << std::endl;  // C++ ostream method

#ifndef sim_LArVoxelID_h
#define sim_LArVoxelID_h

#include <TLorentzVector.h>
#include <TVector3.h>

#include <vector>
#include <functional> // so we can redefine less<> below
#include <iostream>

namespace sim {

  class LArVoxelID 
  {
  public:
    // What is an appropriate storage class for an element of the bin
    // ID?  Probably a regular integer is the only thing we can use; a
    // "short int" can only accept values up to +/-32767, and a large
    // LAr TPC will probably contain more bins on a single axis than
    // that.

    // Construct a voxel identifier given the pre-computed bins.  This
    // one is provided for completeness, but it's not the one I
    // anticipate the clients will use.  Note that it's also serves
    // as the no-argument constructor that ROOT requires for I/O.
    LArVoxelID( const int x = 0, 
		const int y = 0, 
		const int z = 0, 
		const int t = 0 );

    // The constructors that I expect most clients to use: Given the
    // (x,y,z,t) in units of (distance,time) of a particle, compute
    // the voxel identifier.  Note that units are ambiguous there; they
    // have to be consistent between those used by LArVoxelCalculator
    // and the Monte Carlo, but this class does not enforce that
    // consistency.
    explicit LArVoxelID( const TLorentzVector& v );
    LArVoxelID( const double x, 
		const double y, 
		const double z, 
		const double t );

    // Destructor.
    virtual ~LArVoxelID();

  private:
    // Implement the bin ID as a "tuple" of four numbers, one bin for
    // each (x,y,z,t) axis.

    std::vector<int> fbins; // (x,y,z,t) bins.

  public:
    // Accessors for the bin values.  I don't expect these to be used
    // often, but include them for completeness.
    int XBin() const; 
    int YBin() const; 
    int ZBin() const; 
    int TBin() const; 

    // The accessors I expect to be used: The values of the
    // co-ordinates at the bin centers.
    double X() const;
    double Y() const;
    double Z() const;
    double T() const;

    // Alternative accessor: Pretend the LArVoxelID is a vector that
    // takes subscripts, and returns the bin centers of the axis in
    // the square brackets (x=0, y=1, z=2, t=3).
    double operator[]( const int ) const;

    // Conversion function from LArVoxelID to a TLorentzVector.  This
    // lets you do something like:
    //    sim::LArVoxelID larVoxelID = ... // whatever
    //    TLorentzVector v = TLorentzVector( larVoxelID );
    operator TLorentzVector() const;

    // The three-vector version of the above conversion; e.g.:
    //    sim::LArVoxelID larVoxelID = ... // whatever
    //    TLVector3 v = TVector3( larVoxelID );
    operator TVector3() const;

    // The comparison operator.  This a key function, since it
    // establishes the sort order of the voxels in a list.
    bool operator<( const LArVoxelID& ) const;

    // Test for equality.  Handy, but not usually necessary.
    bool operator==( const LArVoxelID& ) const;

    friend std::ostream& operator<< ( std::ostream& output, const LArVoxelID& );

  };

} // sim 

inline int sim::LArVoxelID::XBin() const { return fbins[0]; }
inline int sim::LArVoxelID::YBin() const { return fbins[1]; }
inline int sim::LArVoxelID::ZBin() const { return fbins[2]; }
inline int sim::LArVoxelID::TBin() const { return fbins[3]; }

// A potentially handy definition: At this stage, I'm not sure
// whether I'm going to be keeping a list based on LArVoxelID or on
// LArVoxelID*.  We've already defined operator<(LArVoxelID,LArVoxelID),
// that is, how to compare two LArVoxelID objects; by default that
// also defines std::less<LArVoxelID>, which is what the STL containers
// use for comparisons.

// The following defines std::less<LArVoxelID*>, that is, how to
// compare two LArVoxelID*: by looking at the objects, not at the
// pointer addresses.  The result is that, e.g., a
// map<LArVoxelID*,double> will be sorted in the order I expect.

namespace std {
  template <> 
  class less<sim::LArVoxelID*>
  {
  public:
    bool operator()( const sim::LArVoxelID* lhs, const sim::LArVoxelID* rhs )
    {
      return (*lhs) < (*rhs);
    }
  };
} // std

#endif // sim_LArVoxelID_h
