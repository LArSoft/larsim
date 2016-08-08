////////////////////////////////////////////////////////////////////////
/// $Id: SimChannel.h,v 1.3 2010/03/26 20:08:36 brebel Exp $
///
/// \file  SimChannel.h
///
/// \brief object containing MC truth information necessary for making RawDigits 
/// and doing back tracking
///
/// \author  seligman@nevis.columbia.edu
///
////////////////////////////////////////////////////////////////////////

#ifndef SIM_SIMCHANNEL_H
#define SIM_SIMCHANNEL_H

#include <string>
#include <vector>
#include <map>
#include <stdint.h>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

namespace sim {


  struct TrackIDE{
    int trackID;      ///< Geant4 supplied trackID
    float energyFrac; ///< fraction of hit energy from the particle with this trackID
    float energy;     ///< energy from the particle with this trackID
    TrackIDE() {}
    TrackIDE(int id, float ef, float e) : trackID(id), energyFrac(ef), energy (e) {}
  };

  class IDE{
  public:
    
    IDE();

    //constructor for IDEs applying G4 offset...
    IDE(IDE const&, int);
    
    int    trackID;      ///< Geant4 supplied track ID
    float numElectrons; ///< total number of electrons for this track ID and time
    float energy;       ///< total energy deposited for this track ID and time
    float x;            ///< x position of ionization
    float y;            ///< y position of ionization
    float z;            ///< z position of ionization
  };
  
  class SimChannel
  {
  public:

    // Default constructor
    SimChannel();
    
  private:
    
    raw::ChannelID_t                                    fChannel; ///< electronics channel associated with these sim::Electrons
    std::map< unsigned short, std::vector< sim::IDE > > fTDCIDEs; ///< vector of IDE structs for each TDC with signal


#ifndef __GCCXML__
  public:

    explicit SimChannel(raw::ChannelID_t channel);

    // method to add ionization electrons and energy to this channel
    void AddIonizationElectrons(int trackID,
				unsigned int tdc,
				double numberElectrons,
				double *xyz,
				double energy); 

    
    
    raw::ChannelID_t Channel() const;

    // method to return a collection of IDE structs for all geant4
    // track ids represented between startTDC and endTDC
    std::vector<sim::IDE> TrackIDsAndEnergies(unsigned int startTDC,
					      unsigned int endTDC) const;

    const std::map<unsigned short, std::vector<sim::IDE> >& TDCIDEMap() const;

    // The number of ionization electrons associated with this channel for the 
    // specified TDC.
    double Charge(unsigned int tdc) const;
    double Energy(unsigned int tdc) const;

    // A vector of TrackIDEs for a range of TDCs
    std::vector<sim::TrackIDE> TrackIDEs(unsigned int startTDC,
					 unsigned int endTDC) const;
    
    bool operator<  (const SimChannel& other)     const;
    bool operator== (const SimChannel& other)     const;

    std::pair<int,int> MergeSimChannel(const SimChannel&, int);
    
    //@{
    /**
	  * @brief Dumps the full content of the SimChannel into a stream
	  * @param OSTREAM an ostream-line stream object
	  * @param out the stream to send the information into
	  * @param indent indentation of the lines
	  * @param indent_first indentation for the first line (default: as indent)
	  */
	 template <class OSTREAM>
	 void Dump(OSTREAM& out, std::string indent, std::string first_indent) const;
	 template <class OSTREAM>
	 void Dump(OSTREAM& out, std::string indent = "") const
	   { Dump(out, indent, indent); }
    //@}
    
#endif

  };

} // namespace sim

#ifndef __GCCXML__

inline bool sim::SimChannel::operator<  (const sim::SimChannel& other)                       const { return fChannel < other.Channel(); }
inline bool sim::SimChannel::operator== (const sim::SimChannel& other)                       const { return fChannel == other.Channel(); }
inline const std::map<unsigned short, std::vector<sim::IDE> >& sim::SimChannel::TDCIDEMap() const { return fTDCIDEs; }
inline raw::ChannelID_t sim::SimChannel::Channel()                                          const { return fChannel; }


// -----------------------------------------------------------------------------
// ---  template implementation
// ---
template <class OSTREAM>
void sim::SimChannel::Dump
  (OSTREAM& out, std::string indent, std::string first_indent) const
{
  out << first_indent << "channel #" << Channel() << " read " << fTDCIDEs.size()
    << " TDCs:\n";
  double channel_energy = 0., channel_charge = 0.;
  for (const auto& TDCinfo: fTDCIDEs) {
    unsigned short int tdc = TDCinfo.first;
    out << indent << "  TDC #" << tdc
      << " with " << TDCinfo.second.size() << " IDEs\n";
    double tdc_energy = 0., tdc_charge = 0.;
    for (const sim::IDE& ide: TDCinfo.second) {
      out << indent
        << "    (" << ide.x << ", " << ide.y << ", " << ide.z << ") "
        << ide.numElectrons << " electrons, " << ide.energy << " MeV (trkID="
        << ide.trackID << ")\n";
      tdc_energy += ide.energy;
      tdc_charge += ide.numElectrons;
    } // for IDEs
    out << indent << "    => TDC #" << tdc << " CH #" << Channel()
      << " collected " << tdc_energy << " electrons and " << tdc_energy
      << " MeV\n";
    channel_energy += tdc_energy;
    channel_charge += tdc_charge;
  } // for TDCs
  out << indent << "  => channel #" << Channel() << " collected "
    << channel_charge << " electrons and " << channel_energy << " MeV\n";
} // sim::SimChannel::Dump<>()

#endif

#endif // SIM_SIMCHANNEL_H

////////////////////////////////////////////////////////////////////////
