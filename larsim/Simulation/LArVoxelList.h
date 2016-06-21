////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelList.h
/// \brief Container of LAr voxel information
///
/// \version $Id: LArVoxelList.h,v 1.7 2010/02/15 20:27:06 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// A container for LAr voxel information.  Although there's nothing
/// in the class below that assumes units, the standard for LArSoft is
/// that distances are in cm, and energy is in GeV.
///
/// It acts like a map<LArVoxelID,LArVoxelData>, but with additional
/// features:
///
/// - An Add(LArVoxelList) method allows you to add one LArVoxelList to
///   another.  Voxels with the same ID are added together; otherwise
///   the lists are merged.  Usage:
///      sim::LArVoxelList a,b;
///      a.Add(b);
///   There's even an operator+ so you can write "a += b", but be 
///   careful about wasting memory if you write "a = a+b".  Of
///   course, you can do:
///      sim:LArVoxelList c = a + b;
///
/// - An operator* so you can scale all the voxel energies for
///   calibration:
///      sim::LArVoxelList c;
///      double calib = 1.00000012
///      c *= calib;
///      c = c * calib;  // same as above, but wastes memory
///      sim::LArVoxelList d = c * calib;
///   (Yes, I know we probably won't do calibration in this way,
///   but it's here if you need it.)
///
/// - A method Cut(double) that will remove all voxels with energy
///   less than the argument.
///
/// - Methods ID(int) and Energy(int) for those who are unfamiliar with the
///   concept of "first" and "second" as used with STL maps:
///      sim::LArVoxelList* voxelList = // ...;
///      int numberOfVoxels = voxelList->size();
///      for (int i=0; i<numberOfVoxels; ++i)
///        {
///           sim::LArVoxelID = voxelList->ID(i);
///           double energy = voxelList->Energy(i);
///        }
///   The STL equivalent to the above statements (more efficient):
///      sim::LArVoxelList* voxelList = // ... ;
///      for ( sim::LArVoxelList::const_iterator i = voxelList->begin();
///            i != voxelList->end(); ++i )
///        {
///            sim::LArVoxelID = (*i).first;
///            double energy = (*i).second.Energy();
///        }
///
/// - operator<< method for ROOT display and ease of debugging.

#ifndef LARVOXELLIST_H
#define LARVOXELLIST_H

#include "larsim/Simulation/LArVoxelID.h"
#include "larsim/Simulation/LArVoxelData.h"

#include <map>

namespace sim {

  class LArVoxelList
  {
  public:
    /// Some type definitions to make life easier, and to help "hide"
    /// the implementation details.  (If you're not familiar with STL,
    /// you can ignore these definitions.)
    typedef std::map<LArVoxelID, LArVoxelData>  list_type;
    typedef list_type::key_type                 key_type;
    typedef list_type::mapped_type              mapped_type;
    typedef list_type::value_type               value_type;
    typedef list_type::iterator                 iterator;
    typedef list_type::const_iterator           const_iterator;
    typedef list_type::reverse_iterator         reverse_iterator;
    typedef list_type::const_reverse_iterator   const_reverse_iterator;
    typedef list_type::size_type                size_type;
    typedef list_type::difference_type          difference_type;
    typedef list_type::key_compare              key_compare;
    typedef list_type::allocator_type           allocator_type;

    // Standard constructor and destructor.
    LArVoxelList();
    virtual ~LArVoxelList();

    // Add the energy to the entry with the key; if the key doesn't
    // exist, create it.  Allow the addition both with and without a
    // particle's track ID for LArVoxelData.
    void Add( const key_type& key, const double& energy )                { m_voxelList[key].Add(energy); }
    void Add( const key_type& key, const double& energy, const int& id ) { m_voxelList[key].Add(energy,id); }

    // The arithmetic methods "advertised" in the class description
    // above.
    void Add( const LArVoxelList& );
    LArVoxelList& operator+=( const LArVoxelList& other ) 
    { 
      this->Add(other); 
      return *this; 
    }
    const LArVoxelList operator+(const LArVoxelList& other) const 
    {
      return LArVoxelList(*this) += other;
    }

    LArVoxelList&      operator*=( const double& value );
    const LArVoxelList operator* (const double& value) const 
    {
      return LArVoxelList(*this) *= value;
    }
    // Just in case: define the result of "scalar * LArVoxelList" to be
    // the same as "LArVoxelList * scalar".
    friend const LArVoxelList operator*(const double& value, const LArVoxelList& list);

    // Apply a threshold cut to the voxels in the list, removing all
    // those that fall below the cut.
    void Cut( const double& );

    const key_type& ID( const size_type ) const;
    double Energy( const size_type ) const;

    friend std::ostream& operator<< ( std::ostream& output, const LArVoxelList& );

    // Standard STL methods, to make this class look like an STL map.
    // Again, if you don't know STL, you can just ignore these
    // methods.
    iterator               begin()        { return m_voxelList.begin(); }
    const_iterator         begin()  const { return m_voxelList.begin(); }
    iterator               end()          { return m_voxelList.end(); }
    const_iterator         end()    const { return m_voxelList.end(); }
    reverse_iterator       rbegin()       { return m_voxelList.rbegin(); }
    const_reverse_iterator rbegin() const { return m_voxelList.rbegin(); }
    reverse_iterator       rend()         { return m_voxelList.rend(); }
    const_reverse_iterator rend()   const { return m_voxelList.rend(); }

    size_type size()                const { return m_voxelList.size(); }
    bool      empty()               const { return m_voxelList.empty(); }
    void      swap( LArVoxelList& other ) { m_voxelList.swap( other.m_voxelList ); }
    void      clear()                     { m_voxelList.clear(); }

    iterator       find(const key_type& key)              { return m_voxelList.find(key); }
    const_iterator find(const key_type& key)        const { return m_voxelList.find(key); }
    iterator       upper_bound(const key_type& key)       { return m_voxelList.upper_bound(key); }
    const_iterator upper_bound(const key_type& key) const { return m_voxelList.upper_bound(key); }
    iterator       lower_bound(const key_type& key)       { return m_voxelList.lower_bound(key); }
    const_iterator lower_bound(const key_type& key) const { return m_voxelList.lower_bound(key); }

    mapped_type& operator[](const key_type& key) { return m_voxelList[key]; }
    // My own little addition: operator[] in a const context.
    const mapped_type& operator[]( const key_type& key ) const { return m_voxelList.at(key); }
    mapped_type& at(const key_type& key) { return m_voxelList.at(key); }
    const mapped_type& at(const key_type& key) const { return m_voxelList.at(key); }

    // In addition to operator[], include one insert() method for
    // anyone who wants to include a LArVoxelID / LArVoxelData pair
    // directly.  Note that, as with operator[], there's no check
    // against overwriting an existing item.
    void insert( const key_type& key, const mapped_type& value ) { m_voxelList[key] = value; }

    size_type erase( const key_type& key ) { return m_voxelList.erase(key); }

  private:
    list_type m_voxelList; ///< A sorted list of <LArVoxelID,double> pairs = (voxel ID, energy)

  };
  
} // namespace sim

#endif // LARVOXELLIST_H

