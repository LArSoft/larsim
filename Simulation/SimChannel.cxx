/// $Id: SimChannel.cxx,v 1.3 2010/03/26 20:08:36 brebel Exp $
///
/// \file  Simulation/SimChannel.cxx
///
///
/// \author  seligman@nevis.columbia.edu
///
////////////////////////////////////////////////////////////////////////


#include "Simulation/SimChannel.h"
#include "Simulation/sim.h"
#include "SimpleTypesAndConstants/PhysicalConstants.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace sim{

  //-------------------------------------------------
  IDE::IDE()
    : trackID     (util::kBogusI)
    , numElectrons(util::kBogusD)
    , energy      (util::kBogusD)
    , x           (util::kBogusD)
    , y		  (util::kBogusD)
    , z		  (util::kBogusD)
  {}

  // Default constructor
  //-------------------------------------------------
  SimChannel::SimChannel() 
    : fChannel(0)
  {}

  //-------------------------------------------------
  SimChannel::SimChannel(uint32_t channel)
    : fChannel(channel)
  {}

  
  //-------------------------------------------------
  void SimChannel::AddIonizationElectrons(int trackID,
					  unsigned int tdc,
					  double numberElectrons,
					  double *xyz,
					  double energy)
  {
    // look at the map to see if the current TDC already 
    // exists, if not, add it, if so, just add a new track id to the 
    // vector

    if( fTDCIDEs.count(tdc) > 0 ){      
      // loop over the IDE vector for this tdc and add the electrons 
      // to the entry with the same track id
      std::vector<sim::IDE>::iterator itr = fTDCIDEs[tdc].begin();
      while( itr != fTDCIDEs[tdc].end() ){
	
	if( (*itr).trackID == trackID ){
	  // make a weighted average for the location information
	  double weight       = (*itr).numElectrons + numberElectrons;
	  (*itr).x            = ((*itr).x*(*itr).numElectrons + xyz[0]*numberElectrons)/weight;
	  (*itr).y            = ((*itr).y*(*itr).numElectrons + xyz[1]*numberElectrons)/weight;
	  (*itr).z            = ((*itr).z*(*itr).numElectrons + xyz[2]*numberElectrons)/weight;	  
	  (*itr).numElectrons = weight;
	  (*itr).energy       = (*itr).energy + energy;
	  // found the track id we wanted, so return;
	  return;
	}
	
	itr++;
      }

      // if we never found the track id, then this is the first instance of
      // the track id for this tdc, so add ide to the vector
      sim::IDE ide;
      ide.trackID      = trackID;
      ide.numElectrons = numberElectrons;
      ide.x            = xyz[0];
      ide.y            = xyz[1];
      ide.z            = xyz[2];
      ide.energy       = energy;
      
      fTDCIDEs[tdc].push_back(ide);
    }
    else{
      sim::IDE ide;
      ide.trackID      = trackID;
      ide.numElectrons = numberElectrons;
      ide.x            = xyz[0];
      ide.y            = xyz[1];
      ide.z            = xyz[2];
      ide.energy       = energy;
      
      std::vector<sim::IDE> idelist;
      idelist.push_back(ide);
      fTDCIDEs[tdc] = std::move(idelist);
    }

    return;
  }


  //-------------------------------------------------
  double SimChannel::Charge(unsigned int tdc) const
  {
    double charge = 0.;

    // check to see if this tdc value is in the map
    if( fTDCIDEs.find(tdc) != fTDCIDEs.end() ){
      // loop over the list for this tdc value and add up
      // the total number of electrons
      std::vector<sim::IDE> idelist((*(fTDCIDEs.find(tdc))).second);
      std::vector<sim::IDE>::const_iterator itr = idelist.begin();
      while( itr != idelist.end() ){
	charge += (*itr).numElectrons;
	itr++;
      } // end loop over sim::IDE for this tdc

    } // end if this tdc is represented in the map
  
    return charge;
  }

    //-------------------------------------------------
  double SimChannel::Energy(unsigned int tdc) const
  {
    double energy = 0.;

    // check to see if this tdc value is in the map
    if( fTDCIDEs.find(tdc) != fTDCIDEs.end() ){
      // loop over the list for this tdc value and add up
      // the total number of electrons
      std::vector<sim::IDE> idelist((*(fTDCIDEs.find(tdc))).second);
      std::vector<sim::IDE>::const_iterator itr = idelist.begin();
      while( itr != idelist.end() ){
	energy += (*itr).energy;
	itr++;
      } // end loop over sim::IDE for this tdc

    } // end if this tdc is represented in the map
  
    return energy;
  }

  
  //-----------------------------------------------------------------------
  // the start and end tdc values are assumed to be inclusive
  std::vector<sim::IDE> SimChannel::TrackIDsAndEnergies(unsigned int startTDC,
							unsigned int endTDC) const
  {
    // make a map of track ID values to sim::IDE objects
    std::map<int, sim::IDE> idToIDE;

    std::vector<sim::IDE> ides;

    if(startTDC > endTDC ){
      mf::LogWarning("SimChannel") << "requested tdc range is bogus: "
				   << startTDC << " " << endTDC
				   << " return empty vector";
      return ides;
    }

    std::map<unsigned short, std::vector<sim::IDE> >::const_iterator mitr;
    std::map<unsigned short, std::vector<sim::IDE> >::const_iterator start = fTDCIDEs.lower_bound(startTDC);
    std::map<unsigned short, std::vector<sim::IDE> >::const_iterator end   = fTDCIDEs.upper_bound(endTDC);

    for(mitr = start; mitr != end; mitr++){

      // grab the vector of IDEs for this tdc
      const std::vector<sim::IDE> &idelist = (*mitr).second;
      std::vector<sim::IDE>::const_iterator itr = idelist.begin();
      // now loop over them and add their content to the map
      while( itr != idelist.end() ){
	
	if( idToIDE.find((*itr).trackID) != idToIDE.end() ){
	  double nel1   = idToIDE[(*itr).trackID].numElectrons;
	  double nel2   = (*itr).numElectrons;
	  double en1    = idToIDE[(*itr).trackID].energy;
	  double en2	= (*itr).energy;
	  double energy = en1+en2;
	  double weight = nel1 + nel2;
	  // make a weighted average for the location information
	  idToIDE[(*itr).trackID].x            = ((*itr).x*nel2 + idToIDE[(*itr).trackID].x*nel1)/weight;
	  idToIDE[(*itr).trackID].y            = ((*itr).y*nel2 + idToIDE[(*itr).trackID].y*nel1)/weight;
	  idToIDE[(*itr).trackID].z            = ((*itr).z*nel2 + idToIDE[(*itr).trackID].z*nel1)/weight;	  
	  idToIDE[(*itr).trackID].numElectrons = weight;
	  idToIDE[(*itr).trackID].energy = energy;
	} // end if the track id for this one is found
	else{
	  sim::IDE temp(*itr);
	  idToIDE[(*itr).trackID] = temp;
	}

	itr++;
      } // end loop over vector
    } // end loop over tdc values

    // now fill the vector with the ides from the map
    for(std::map<int, sim::IDE>::iterator itr = idToIDE.begin(); itr != idToIDE.end(); itr++){
      ides.push_back((*itr).second);
    }

    return ides;
  }

  //-----------------------------------------------------------------------
  // the start and end tdc values are assumed to be inclusive
  std::vector<sim::TrackIDE>  SimChannel::TrackIDEs(unsigned int startTDC,
						      unsigned int endTDC) const
  {

    std::vector<sim::TrackIDE> trackIDEs;

    if(startTDC > endTDC ){
      mf::LogWarning("SimChannel::TrackIDEs") << "requested tdc range is bogus: "
					      << startTDC << " " << endTDC
					      << " return empty vector";
      return trackIDEs;
    }

    double totalE = 0.;
    const std::vector<sim::IDE> ides = TrackIDsAndEnergies(startTDC,endTDC);    
    for (auto const& ide : ides)
      totalE += ide.energy;

    // protect against a divide by zero below
    if(totalE < 1.e-5) totalE = 1.;
    
    // loop over the entries in the map and fill the input vectors    
    for (auto const& ide : ides){      
      if(ide.trackID == sim::NoParticleId) continue;
      trackIDEs.emplace_back(ide.trackID,ide.energy/totalE,ide.energy); 
    }


    return trackIDEs;
  }

}
