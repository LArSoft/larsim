// SimPhotons, OnePhoton and SimPhotonsCollection implementation.
//
// These objects are primarily storage containers, so not much
// function code is actually required.
//
// See comments at head of SimPhotons.h for more details.
//
// Ben Jones, MIT, 06/04/2010
//
/// \version $Id: ParticleHistory.cxx,v 1.1 2010/04/29 15:38:01 seligman Exp $


#include "Simulation/SimPhotons.h"
#include <iostream>

namespace sim
{

  //----------------------------------------------------------------------------
  SimPhotonsCollection::SimPhotonsCollection()
  {
  }

  //----------------------------------------------------------------------------
  /*
    SimPhotons SimPhotonsCollection::GetHit(int key)
  {
    if( !((*this)[key]) )
      (*this)[key] = new SimPhotons();
     return (*this)[key];
  }
  */ 
  //----------------------------------------------------------------------------
  /*
  SimPhotonsCollection & SimPhotonsCollection::operator+=(const SimPhotonsCollection &rhs)
  {
    for(SimPhotonsCollection::const_iterator it = rhs.begin(); it!=rhs.end(); it++){
      GetHit(it->first).operator+=(it->second);
    }
    SetSDName("CompositeHitCollection");
    return *this;
  }
  */
  //----------------------------------------------------------------------------
  /*
  const SimPhotonsCollection SimPhotonsCollection::operator+(const SimPhotonsCollection &rhs) const
  {
    return SimPhotonsCollection(*this)+=rhs;
  }
  */
  //----------------------------------------------------------------------------
  OnePhoton::OnePhoton()
  {
  }

  //----------------------------------------------------------------------------
  SimPhotonsLite::SimPhotonsLite()
  {
  }

  /// \todo: Remove this constructor when LBNE makes the next round of production
  ///        MC files - after 11 September 2013 brebel
  //----------------------------------------------------------------------------
  LBNE10ktPhotons::LBNE10ktPhotons()
  {
  }

  //----------------------------------------------------------------------------
  SimPhotons::SimPhotons()
  {
  }

  //----------------------------------------------------------------------------
  SimPhotons & SimPhotons::operator+=(const SimPhotons &rhs)
  {
    for(SimPhotons::const_iterator it = rhs.begin(); it!=rhs.end(); it++){
      push_back(*it);
    }
    return *this;
  }

  //----------------------------------------------------------------------------
  const SimPhotons SimPhotons::operator+(const SimPhotons &rhs) const
  {
    return SimPhotons(*this)+=rhs;
  }

}
