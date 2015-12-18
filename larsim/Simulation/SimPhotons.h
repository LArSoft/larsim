////////////////////////////////////////////////////////////////////////
/// \file  SimPhotons.h
/// \brief contains objects relating to OpDet hits
///
/// \version $Id: ParticleList.h,v 1.13 2010/05/13 16:12:20 seligman Exp $
/// \author  Ben Jones
////////////////////////////////////////////////////////////////////////
// This file contains the definitions of the classes which
// are stored in the event representing OpDet hits.
//
// A OpDet Hit stores data for each photon which steps inside the OpDet 
// volume.  Currently the quantities stored are 4 potition, 4 momentum
// and TrackID.  A SimPhotonsCollection is a set of SimPhotonss, one per OpDet
// in the collection. 
// 
// The SimPhotons is filled in by the OpDetSensitiveDetector class in LArG4 
// and will be used to generate the OpDet response later in the simulation
// chain.
//
// OnePhoton, SimPhotons and SimPhotonsCollection are all persistent under
// ROOT I/O.  For compilation to succeed, the relevant pragma lines
// must be present in LinkDef.h.
//
// The current implementation resembles that of an STL container in 
// some respects but needs more work before it is polished product.
//
// Ben Jones, MIT, 06/04/2010
//

#ifndef SimPhotons_h
#define SimPhotons_h

#include "TLorentzVector.h"

#include <map>


namespace sim
{

  // This structure contains all the information per photon
  // which entered the sensitive OpDet volume.

  class OnePhoton 
  {
  public:
    OnePhoton();

    bool           SetInSD;
    TVector3       InitialPosition;
    TVector3       FinalLocalPosition; // in cm
    float          Time;
    float          Energy;
  };
  
  class SimPhotonsLite
  {
    public:
      SimPhotonsLite();
      SimPhotonsLite(int chan)
	: OpChannel(chan)
      {}
      
      int   OpChannel;
      std::map<int, int> DetectedPhotons;

      SimPhotonsLite& operator+=(const SimPhotonsLite &rhs);
      const SimPhotonsLite operator+(const SimPhotonsLite &rhs) const;

      bool operator==(const SimPhotonsLite &other) const;
  };

  /// \todo: Remove this class when DUNE makes the next round of production
  ///        MC files - after 11 September 2013 brebel
  class DUNE10ktPhotons
  {
    public:
      DUNE10ktPhotons();
      int   OpChannel;
      std::map<int, int> DetectedPhotons;
  };

  // Define a OpDet Hit as a list of OpDet photons which were
  // recorded in the OpDet volume.

  class SimPhotons : public std::vector<OnePhoton> 
    {
    public:
      SimPhotons();
      SimPhotons(int chan)
	: fOpChannel(chan)
      {}

      int  fOpChannel;  /// volume number for the OpDet

#ifndef __GCCXML__
    public:

      typedef std::vector<OnePhoton>             list_type;
      typedef list_type::value_type              value_type;
      typedef list_type::iterator                iterator;
      typedef list_type::const_iterator          const_iterator;
      typedef list_type::reverse_iterator        reverse_iterator;
      typedef list_type::const_reverse_iterator  const_reverse_iterator;
      typedef list_type::size_type               size_type;
      typedef list_type::difference_type         difference_type;

      // define addition operators for combining hits
      //   (add all photons to one vector)
      SimPhotons& operator+=(const SimPhotons &rhs);
      const SimPhotons operator+(const SimPhotons &rhs) const;

      bool operator== (const SimPhotons &other) const;
      
      int       OpChannel() const;
      void      SetChannel(int ch);

#endif
      
    };
 


  // The OpDet Hit collection is the set of all OpDet Hits indexed
  // by OpDet ID 

  class SimPhotonsCollection : public std::map<int, SimPhotons>{
  public:

    SimPhotonsCollection();

  private:
    std::string fTheSDName;
    
#ifndef __GCCXML__
    
  public:
    typedef std::map<int,SimPhotons>           list_type;
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
    
    //SimPhotons&  GetHit(int);
    //SimPhotons  GetHit(int);
    
    // define addition operators for combining hit collections
    //   (add each hit in the collection)
    //SimPhotonsCollection& operator+=(const SimPhotonsCollection &rhs);
    //const SimPhotonsCollection operator+(const SimPhotonsCollection &rhs) const; 
    
  public:
    void SetSDName(std::string TheSDName);
    std::string GetSDName();
    
#endif
  };
  
}

#ifndef __GCCXML__

inline int         sim::SimPhotons::OpChannel()       const                     { return fOpChannel;      }
inline void        sim::SimPhotons::SetChannel(int ch)                          { fOpChannel = ch;        }
inline std::string sim::SimPhotonsCollection::GetSDName()                       { return fTheSDName;      }
inline void        sim::SimPhotonsCollection::SetSDName(std::string TheSDName)  { fTheSDName = TheSDName; }

inline bool sim::SimPhotons::operator==(const sim::SimPhotons& other) const          { return fOpChannel == other.OpChannel(); }
inline bool sim::SimPhotonsLite::operator==(const sim::SimPhotonsLite& other) const  { return OpChannel == other.OpChannel; }

#endif

#endif
