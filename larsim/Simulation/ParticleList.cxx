////////////////////////////////////////////////////////////////////////
/// \file  ParticleList.cxx
/// \brief Particle list in DetSim contains Monte Carlo particle information.
///
/// \version $Id: ParticleList.cxx,v 1.10 2010/05/13 16:12:20 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// Although there's nothing in the following class that assumes
/// units, the standard for LArSoft is that distances are in cm, and
/// energies are in GeV.

#include "larsim/Simulation/ParticleList.h"
#include "larsim/Simulation/EveIdCalculator.h"

#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <TLorentzVector.h>

#include <set>
// #include <iterator>
//#include <cmath>
// #include <memory>

namespace sim {

  //----------------------------------------------------------------------------
  // Constructor.
  ParticleList::ParticleList() 
  {
  }  

  //----------------------------------------------------------------------------
  // Destructor 
  ParticleList::~ParticleList()
  {
    this->clear();
  }

  //----------------------------------------------------------------------------
  // Copy constructor.  Note that since this class inherits from
  // TObject, we have to copy its information explicitly.
  ParticleList::ParticleList( const ParticleList& rhs ) 
  {
    // Clear any contents that we already possess.
    this->clear();

    // Copy each entry in the other ParticleList.
    for ( const_iterator entry = rhs.m_particleList.begin();
	  entry != rhs.m_particleList.end(); ++entry ){
      const simb::MCParticle* original = (*entry).second;
      simb::MCParticle* copy = new simb::MCParticle( *original );
      this->insert( copy );
    }
  }

  //----------------------------------------------------------------------------
  // Assignment constructor.
  ParticleList& ParticleList::operator=( const ParticleList& rhs )
  {
    // Usual test for self-assignment.
    if ( this == &rhs ) return *this;

    // Clear any contents that we already possess.
    this->clear();

    // Copy each entry in the other ParticleList.
    for ( const_iterator entry = rhs.m_particleList.begin();
	  entry != rhs.m_particleList.end(); ++entry ){
      const simb::MCParticle* original = (*entry).second;
      simb::MCParticle* copy = new simb::MCParticle( *original );
      this->insert( copy );
    }
    
    return *this;
  }

  //----------------------------------------------------------------------------
  // Apply an energy cut to the particles.
  void ParticleList::Cut( const double& cut )
  {
    // The safest way to do this is to create a list of track IDs that
    // fail the cut, then delete those IDs.

    // Define a list of IDs.
    typedef std::set< key_type > keyList_type;
    keyList_type keyList;

    // Add each ID that fails the cut to the list.
    for ( const_iterator i = m_particleList.begin(); i != m_particleList.end(); ++i ){
      const simb::MCParticle* particle = (*i).second;
      Double_t totalInitialEnergy = particle->E();
      if ( totalInitialEnergy < cut ) { 
	keyList.insert( (*i).first ); 
      }
    }
    
    // Go through the list, deleting the particles that are on the list.
    for ( keyList_type::const_iterator i = keyList.begin(); i != keyList.end(); ++i ){
      this->erase( *i );
    }
  }


  //----------------------------------------------------------------------------
  const ParticleList::key_type& ParticleList::TrackId( const size_type index ) const
  {
    const_iterator i = m_particleList.begin();
    std::advance(i,index);
    return (*i).first;
  }
  //----------------------------------------------------------------------------
  ParticleList::mapped_type ParticleList::Particle( const size_type index ) const
  {
    const_iterator i = m_particleList.begin();
    std::advance(i,index);
    return (*i).second;
  }

  //----------------------------------------------------------------------------
  ParticleList::mapped_type ParticleList::Particle( const size_type index )
  {
    iterator i = m_particleList.begin();
    std::advance(i,index);
    return (*i).second;
  }

  //----------------------------------------------------------------------------
  bool ParticleList::IsPrimary( int trackID ) const
  {
    return m_primaries.find( trackID )  !=  m_primaries.end();
  }

  //----------------------------------------------------------------------------
  int ParticleList::NumberOfPrimaries() const
  {
    return m_primaries.size();
  }

  //----------------------------------------------------------------------------
  const simb::MCParticle* ParticleList::Primary( const int index ) const
  {
    // Advance "index" entries from the beginning of the primary list.
    primaries_const_iterator primary = m_primaries.begin();
    std::advance( primary, index );

    // Get the track ID from that entry in the list.
    int trackID = *primary;

    // Find the entry in the particle list with that track ID.
    const_iterator entry = m_particleList.find(trackID);

    // Return the Particle object in that entry.
    return (*entry).second;
  }

  //----------------------------------------------------------------------------
  simb::MCParticle* ParticleList::Primary( const int index )
  {
    // Advance "index" entries from the beginning of the primary list.
    primaries_const_iterator primary = m_primaries.begin();
    std::advance( primary, index );

    // Get the track ID from that entry in the list.
    int trackID = *primary;

    // Find the entry in the particle list with that track ID.
    iterator entry = m_particleList.find(trackID);

    // Return the Particle object in that entry.
    return (*entry).second;
  }


  //----------------------------------------------------------------------------
  std::vector<const simb::MCParticle*> ParticleList::GetPrimaries() const
  {
    std::vector<const simb::MCParticle*> primaries;
    primaries.reserve(m_primaries.size());
    // for each particle, check if its track ID is in the primaries list
    // iPartPair is std::pair<const int, simb::MCParticle*>
    for (auto& iPartPair: m_particleList) {
      if (m_primaries.count(iPartPair.first))
        primaries.push_back(iPartPair.second);
    } // for
    if (primaries.size() != m_primaries.size()) {
      throw cet::exception("ParticleList")
        << "sim::ParticleList::GetPrimaries() collected " << primaries.size()
        << " primaries, not " << m_primaries.size() << " as expected\n";
    }
    return primaries;
  } // ParticleList::GetPrimaries()
  
  
  //----------------------------------------------------------------------------
//   void ParticleList::Add( const ParticleList& other )
//   {
//     int offset = 0;
//     if ( ! m_particleList.empty() ){
//       // Get the ID number of the last track in our list of particles.
//       const_iterator last = m_particleList.end();
//       --last;
//       int lastTrackID = (*last).first;
      
//       // Compute an offset so that there will be no overlaps in
//       // track numbers.
//       offset = lastTrackID + 1;
      
//       // The first track ID number in the other list is probably 1,
//       // or perhaps 0.  But there's a chance that it might be a
//       // negative number, if the user manually adds a negative track
//       // ID as some indicator.  Let's allow for that.
//       int firstOtherTrackID = 0;
//       if ( ! other.m_particleList.empty() ){
// 	// Get the first track number of the other list.
// 	firstOtherTrackID = (*(m_particleList.begin())).first;
	
// 	if ( firstOtherTrackID < 0 ){
// 	  offset -= firstOtherTrackID;
// 	}
//       }
//     }

//     // Create a new particle list from "other", with non-overlapping
//     // track IDs with respect to our list.
//     ParticleList adjusted = other + offset;
    
//     // Merge the two particle lists.
//     m_particleList.insert( adjusted.m_particleList.begin(), adjusted.m_particleList.end() );
//     m_primaries.insert( adjusted.m_primaries.begin(), adjusted.m_primaries.end() );
//   }

  //----------------------------------------------------------------------------
//   ParticleList ParticleList::Add( const int& offset ) const
//   {
//     // Start with a fresh ParticleList, the destination of the
//     // particles with adjusted track numbers.
//     ParticleList result;

//     // For each particle in our list:
//     for ( const_iterator i = m_particleList.begin(); i != m_particleList.end(); ++i ){
//       const simb::MCParticle* particle = (*i).second;
      
//       // Create a new particle with an adjusted track ID.
//       simb::MCParticle* adjusted = new simb::MCParticle( particle->TrackId() + offset,
// 						   particle->PdgCode(),
// 						   particle->Process(),
// 						   particle->Mother(),
// 						   particle->Mass() );
      
//       adjusted->SetPolarization( particle->Polarization() );
      
//       // Copy all the daughters, adjusting the track ID.
//       for ( int d = 0; d < particle->NumberDaughters(); ++d ){
// 	int daughterID = particle->Daughter(d);
// 	adjusted->AddDaughter( daughterID + offset );
//       }
      
//       // Copy the trajectory points.
//       for ( size_t t = 0; t < particle->NumberTrajectoryPoints(); ++t ){
// 	adjusted->AddTrajectoryPoint( particle->Position(t), particle->Momentum(t) );
//       }
      
//       // Add the adjusted particle to the destination particle list.
//       // This will also adjust the destination's list of primary
//       // particles, if needed.
//       result.insert( adjusted );
//     }
    
//     return result;
//   }
  
  //----------------------------------------------------------------------------
  // Just in case: define the result of "scalar * ParticleList" to be
  // the same as "ParticleList * scalar".
//   ParticleList operator+(const int& value, const ParticleList& list) 
//   {
//     return list + value;
//   }


  // This is the main "insertion" method for the ParticleList
  // pseudo-array pseudo-map.  It does the following:
  //  - Add the Particle to the list; if the track ID is already in the
  //    list, throw an exception.
  //  - If it's a primary particle, add it to the list of primaries.
  void ParticleList::insert( simb::MCParticle* particle ) 
  { 
    int trackID = particle->TrackId();
    iterator insertion = m_particleList.lower_bound( trackID );
    if ( insertion == m_particleList.end() ){
      // The best "hint" we can give is that the particle will go at
      // the end of the list.
      m_particleList.insert( insertion, value_type( trackID, particle ) );
    }
    else if ( (*insertion).first == trackID ){
      throw cet::exception("ParticleList") << "sim::ParticleList::insert - ERROR - "
				   << "track ID=" << trackID 
				   << " is already in the list\n";
    }
    else{
      // It turns out that the best hint we can give is one more
      // than the result of lower_bound.
      m_particleList.insert( ++insertion, value_type( trackID, particle ) );
    }

    // If this is a primary particle, add it to the list.  use 
    // rimary as the process string to look for as the P may or may not
    // be capitalized
    if ( particle->Process().find("rimary") != std::string::npos )
      m_primaries.insert( trackID );
  }

  //----------------------------------------------------------------------------
  void ParticleList::clear()
  {
    for ( iterator i = m_particleList.begin(); i != m_particleList.end(); ++i ){
      delete (*i).second;
    }

    m_particleList.clear();
  }

  //----------------------------------------------------------------------------
  // An erase that includes the deletion of the associated Particle*.
  ParticleList::size_type ParticleList::erase( const key_type& key )
  {
    iterator entry = m_particleList.find( abs(key) );
    delete (*entry).second;
    return m_particleList.erase( abs(key) );
  }


  //----------------------------------------------------------------------------
  std::ostream& operator<< ( std::ostream& output, const ParticleList& list )
  {
    // Determine a field width for the particle number.
    ParticleList::size_type numberOfParticles = list.size();
    int numberOfDigits = (int) std::log10( (double) numberOfParticles ) + 1;

    // A simple header.
    output.width( numberOfDigits );
    output << "#" << ": < ID, particle >" << std::endl; 

    // Write each particle on a separate line.
    ParticleList::size_type nParticle = 0;
    for ( ParticleList::const_iterator particle = list.begin(); 
	  particle != list.end(); ++particle, ++nParticle ){
      output.width( numberOfDigits );
      output << nParticle << ": " 
	     << "<" << (*particle).first 
	     << "," << *((*particle).second) 
	     << ">" << std::endl;
    }

    return output;
  }

  //----------------------------------------------------------------------------
  // operator[] in a non-const context.
  ParticleList::mapped_type ParticleList::operator[]( const key_type& key)
  {
    const_iterator i = m_particleList.find(abs(key));
    if( i == m_particleList.end() ) 
      throw cet::exception("ParticleList") << "i == m_particleList.end()\n";

    return (*i).second;
  }

  // operator[] in a const context; not usual for maps, but I find it useful.
  ParticleList::mapped_type ParticleList::operator[]( const key_type& key) const
  {
    const_iterator i = m_particleList.find(abs(key));
    if( i == m_particleList.end() ) 
      throw cet::exception("ParticleList") << "i == m_particleList.end()\n";

    return (*i).second;
  }



  //----------------------------------------------------------------------------
  // Static variable for eve ID calculator. We're using an unique_ptr,
  // so when this object is eventually deleted (at the end of the job)
  // it will delete the underlying pointer.
  static std::unique_ptr<EveIdCalculator> eveIdCalculator;

  //----------------------------------------------------------------------------
  // The eve ID calculation.
  int ParticleList::EveId( const int trackID ) const
  {
    // If the eve ID calculator has never been initialized, use the
    // default method.
    if ( eveIdCalculator.get() == 0 ){
      AdoptEveIdCalculator( new EveIdCalculator );
    }

    // If the eve ID calculator has changed, or we're looking at a
    // different ParticleList, initialize the calculator.
    static EveIdCalculator* saveEveIdCalculator = 0;
    if ( saveEveIdCalculator != eveIdCalculator.get() ) {
      saveEveIdCalculator = eveIdCalculator.get();
      eveIdCalculator->Init( this );
    }
    if ( eveIdCalculator->ParticleList() != this ){
      eveIdCalculator->Init( this );
    }
    
    // After the "bookkeeping" tests, here's where we actually do the
    // calculation.
    return eveIdCalculator->CalculateEveId( trackID );
  }

  //----------------------------------------------------------------------------
  // Save a new eve ID calculation method.
  void ParticleList::AdoptEveIdCalculator( EveIdCalculator* calc )
  {
    eveIdCalculator.reset(calc);
  }

} // namespace sim
