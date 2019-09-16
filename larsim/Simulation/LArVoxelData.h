////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelData.h
/// \brief Encapsulates the information we want store for a voxel.
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// LArVoxelList associates a voxel ID with voxel data. LArVoxelID
/// describes the ID (position); this class describes the data.
///
/// In particular, this class stores the association between particle
/// tracks and the energy deposited within a voxel.
///
/// The key item of information stored for a voxel is the energy.
/// Another, less important item is the association of a particular
/// particle track to an amount of energy stored in the voxel.  If all
/// you want is the energy in a voxel, use LArVoxelData::Energy() and
/// ignore the rest of this discussion. If you want to understand the
/// source of that energy, you have to dig deeper:
///
/// LArVoxelData::Energy[trackID] returns the amount of energy
/// deposited by the particle with the given track ID (note the use of
/// square brackets).
///
/// LArVoxelData::NumberParticles() returns the number of individual
/// particles whose energy deposits are recorded for this voxel.
///
/// LArVoxelData::Energy(i) returns the amount of energy deposited by
/// the "i-th" particle (note the use of parenthesis).
///
/// LArVoxelData::TrackID(i) returns the track ID of the "i-th"
/// particle.
///
/// LArVoxelData::UnassignedEnergy() returns the amount of energy in
/// the voxel that's not assigned to any particle. There might be
/// "unassigned" energy for one of two reasons:
///
/// 1) In the Monte Carlo, particles can be cut from the ParticleList
///    because they fall below the energy cut. If that happens, the
///    particle is also removed from the LArVoxelData lists, and its
///    energy placed in "unassigned" energy.
///
///    For example, a voxel might contain the sum of many low-energy
///    electrons. It may be that none of the individual electron
///    tracks would be written to the output, but their sum might be
///    enough for the voxel to be written to the output file.
///
/// 2) If any form of "voxel arithmetic" is performed (e.g,. adding
///    two LArVoxelLists together), it becomes too difficult to
///    maintain the particle<->energy relationship because the track
///    IDs change. In that case, all the voxel energies are moved to
///    "unassigned" energy.
///
/// LArVoxelData::AssignedEnergy() returns, in effect,
///    Energy() - UnassignedEnergy()
///

#ifndef Simulation_LArVoxelData_h_
#define Simulation_LArVoxelData_h_

#include "lardataobj/Utilities/VectorMap.h"
#include "lardataobj/Utilities/SumSecondFunction.h"
#include "larsim/Simulation/LArVoxelID.h"

#include <iosfwd>
#include <numeric>

namespace sim {

  class LArVoxelData
  {
  public:
    // Some type definitions to make life easier, and to help "hide"
    // the implementation details.  (If you're not familiar with STL,
    // you can ignore these definitions.)
    typedef util::VectorMap<int, double>       list_type;
    typedef list_type::key_type                key_type;
    typedef list_type::mapped_type             mapped_type;
    typedef list_type::value_type              value_type;
    typedef list_type::iterator                iterator;
    typedef list_type::const_iterator          const_iterator;
    typedef list_type::reverse_iterator        reverse_iterator;
    typedef list_type::const_reverse_iterator  const_reverse_iterator;
    typedef list_type::size_type               size_type;
    typedef list_type::difference_type         difference_type;
    typedef list_type::key_compare             key_compare;
    typedef list_type::allocator_type          allocator_type;

    // Standard constructor and destructor.
    LArVoxelData();
    virtual ~LArVoxelData();

  private:
    // This is the sum of all the voxel energy that is not assigned
    // to a particular particle track.
    mapped_type fenergy;     // Energy not assigned to a particular track in a voxel.

    // If we're able to maintain a track<->energy relationship for
    // this voxel, this map contains the amount of energy deposited
    // in this voxel for the given tracks.
    list_type ftrackEnergy;  // Energy assigned to individual particle tracks in a voxel.

    sim::LArVoxelID fVoxelID; //id for the voxel represented by these data

  public:
    // The energy routines described above.  (std::accumulate is
    // defined in <numeric>, and is a standard STL algorithm.)
    mapped_type AssignedEnergy() const;
    mapped_type UnassignedEnergy() const;
    mapped_type Energy() const;

    size_type NumberParticles() const;

    const key_type&    TrackID( const size_type ) const;
    const mapped_type& Energy ( const size_type ) const;

    // Add the energy to the entry with the key; if the key doesn't
    // exist, create it.
    void Add( const mapped_type& energy, const key_type& trackID );

    // If there's no key, it must be "unassigned" energy.
    void Add( const mapped_type& energy );

    sim::LArVoxelID VoxelID()  const;
    void            SetVoxelID(sim::LArVoxelID voxID);

    // Arithmetic methods to support the arithmetic that can be
    // performed by LArVoxelList.
    void Add( const LArVoxelData& );
    LArVoxelData& operator+=( const LArVoxelData& other );
    const LArVoxelData operator+(const LArVoxelData& other) const;

    LArVoxelData& operator*=( const double& value );
    const LArVoxelData operator*(const double& value) const;

    // Just in case: define the result of "scalar * LArVoxelData" to be
    // the same as "LArVoxelData * scalar".
    friend const LArVoxelData operator*(const double& value, const LArVoxelData& list);

    // RemoveTrack any particle with this track ID and move its energy to
    // the "unassigned" energy; it returns the number of entries
    // removed.  This method is inlined because we need to do this
    // fast; it's probably being called in a loop over all the
    // members of a LArVoxelList.
    size_type RemoveTrack( const int& track );

    // Remove all particles and move their energies to "unassigned."
    void RemoveAllTracks();

    friend std::ostream& operator<< ( std::ostream& output, const LArVoxelData& );

    // Standard STL methods, to make this class look like an STL map.
    // Again, if you don't know STL, you can just ignore these
    // methods. Remember, the "map" portion of this class doesn't
    // always tell the whole story; you also need to look at the
    // "unasigned" energy separately.
    iterator               begin();
    const_iterator         begin()  const;
    iterator               end();
    const_iterator         end()    const;
    reverse_iterator       rbegin();
    const_reverse_iterator rbegin() const;
    reverse_iterator       rend();
    const_reverse_iterator rend()   const;

    size_type size()                const;
    bool      empty()               const;
    void      swap( LArVoxelData& other );
    void      clear();

    iterator       find(const key_type& key);
    const_iterator find(const key_type& key)        const;
    iterator       upper_bound(const key_type& key);
    const_iterator upper_bound(const key_type& key) const;
    iterator       lower_bound(const key_type& key);
    const_iterator lower_bound(const key_type& key) const;
    size_type      erase( const key_type& key );

    mapped_type&       operator[](const key_type& key);
    // My own little addition: operator[] in a const context.
    const mapped_type& operator[]( const key_type& key) const;
    mapped_type&       at(const key_type& key);
    const mapped_type& at(const key_type& key)          const;

    // In addition to operator[], include one insert() method.  Note
    // that, as with operator[], there's no check against overwriting
    // an existing item.
    void insert( const key_type& key, const mapped_type& value );

  };

} // namespace sim

inline void                  sim::LArVoxelData::SetVoxelID(sim::LArVoxelID voxID) { fVoxelID = voxID; }
inline sim::LArVoxelID                sim::LArVoxelData::VoxelID()          const { return fVoxelID;  }
inline sim::LArVoxelData::size_type   sim::LArVoxelData::NumberParticles()  const { return size();    }
inline sim::LArVoxelData::mapped_type sim::LArVoxelData::UnassignedEnergy() const { return fenergy;   }
inline sim::LArVoxelData::mapped_type sim::LArVoxelData::AssignedEnergy()   const
{
  return std::accumulate( ftrackEnergy.begin(), ftrackEnergy.end(), 0.0,
			  util::SumSecondFunction<key_type,mapped_type>() );
}
inline sim::LArVoxelData::mapped_type sim::LArVoxelData::Energy()           const
{
  return std::accumulate( ftrackEnergy.begin(), ftrackEnergy.end(), fenergy,
			  util::SumSecondFunction<key_type,mapped_type>() );
}
inline sim::LArVoxelData::size_type   sim::LArVoxelData::RemoveTrack( const int& track )
{
  iterator search = ftrackEnergy.find( track );
  if ( search != ftrackEnergy.end() )
    {
      fenergy += (*search).second;
      ftrackEnergy.erase( search );
      return 1;
    }
  else return 0;
}
inline void sim::LArVoxelData::RemoveAllTracks()
{
  fenergy = this->Energy();
  ftrackEnergy.clear();
}
inline void sim::LArVoxelData::Add( const sim::LArVoxelData::mapped_type& energy,
				    const sim::LArVoxelData::key_type& trackID )
{ ftrackEnergy[trackID] += energy; }
inline void sim::LArVoxelData::Add( const sim::LArVoxelData::mapped_type& energy)
{ fenergy += energy; }
inline sim::LArVoxelData& sim::LArVoxelData::operator+=( const sim::LArVoxelData& other)
{
  this->Add(other);
  return *this;
}
inline const sim::LArVoxelData sim::LArVoxelData::operator+(const sim::LArVoxelData& other) const
{
  return LArVoxelData(*this) += other;
}
inline const sim::LArVoxelData sim::LArVoxelData::operator*(const double& value) const
{
  return LArVoxelData(*this) *= value;
}
inline sim::LArVoxelData::iterator               sim::LArVoxelData::begin()        { return ftrackEnergy.begin();  }
inline sim::LArVoxelData::const_iterator         sim::LArVoxelData::begin()  const { return ftrackEnergy.begin();  }
inline sim::LArVoxelData::iterator               sim::LArVoxelData::end()          { return ftrackEnergy.end();    }
inline sim::LArVoxelData::const_iterator         sim::LArVoxelData::end()    const { return ftrackEnergy.end();    }
inline sim::LArVoxelData::reverse_iterator       sim::LArVoxelData::rbegin()       { return ftrackEnergy.rbegin(); }
inline sim::LArVoxelData::const_reverse_iterator sim::LArVoxelData::rbegin() const { return ftrackEnergy.rbegin(); }
inline sim::LArVoxelData::reverse_iterator       sim::LArVoxelData::rend()         { return ftrackEnergy.rend();   }
inline sim::LArVoxelData::const_reverse_iterator sim::LArVoxelData::rend()   const { return ftrackEnergy.rend();   }

inline sim::LArVoxelData::size_type              sim::LArVoxelData::size()   const { return ftrackEnergy.size();   }
inline bool                                      sim::LArVoxelData::empty()  const { return ftrackEnergy.empty();  }
inline void sim::LArVoxelData::swap( LArVoxelData& other )
{
  ftrackEnergy.swap( other.ftrackEnergy );
  double temp = fenergy;
  fenergy = other.fenergy;
  other.fenergy = temp;
}
inline void sim::LArVoxelData::clear() { fenergy = 0.; ftrackEnergy.clear(); }
inline sim::LArVoxelData::iterator       sim::LArVoxelData::find(const sim::LArVoxelData::key_type& key)
{ return ftrackEnergy.find(key);        }
inline sim::LArVoxelData::const_iterator sim::LArVoxelData::find(const sim::LArVoxelData::key_type& key) const
{ return ftrackEnergy.find(key);        }
inline sim::LArVoxelData::iterator       sim::LArVoxelData::upper_bound(const sim::LArVoxelData::key_type& key)
{ return ftrackEnergy.upper_bound(key); }
inline sim::LArVoxelData::const_iterator sim::LArVoxelData::upper_bound(const sim::LArVoxelData::key_type& key) const
{ return ftrackEnergy.upper_bound(key); }
inline sim::LArVoxelData::iterator       sim::LArVoxelData::lower_bound(const sim::LArVoxelData::key_type& key)
{ return ftrackEnergy.lower_bound(key); }
inline sim::LArVoxelData::const_iterator sim::LArVoxelData::lower_bound(const sim::LArVoxelData::key_type& key) const
{ return ftrackEnergy.lower_bound(key); }
inline sim::LArVoxelData::size_type      sim::LArVoxelData::erase( const sim::LArVoxelData::key_type& key )
{ return this->RemoveTrack(key);        }

inline sim::LArVoxelData::mapped_type&   sim::LArVoxelData::operator[](const sim::LArVoxelData::key_type& key)
{ return ftrackEnergy[key];    }
inline const sim::LArVoxelData::mapped_type& sim::LArVoxelData::operator[]( const sim::LArVoxelData::key_type& key) const
{ return ftrackEnergy.at(key); }
inline sim::LArVoxelData::mapped_type&       sim::LArVoxelData::at(const sim::LArVoxelData::key_type& key)
{ return ftrackEnergy.at(key); }
inline const sim::LArVoxelData::mapped_type& sim::LArVoxelData::at(const sim::LArVoxelData::key_type& key) const
{ return ftrackEnergy.at(key); }
inline void sim::LArVoxelData::insert( const sim::LArVoxelData::key_type& key, const sim::LArVoxelData::mapped_type& value )
{ ftrackEnergy[key] = value; }

#endif // Simulation_LArVoxelData_h_
