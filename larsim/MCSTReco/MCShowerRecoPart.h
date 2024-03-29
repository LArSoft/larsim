#ifndef MCSHOWERRECOPART_H
#define MCSHOWERRECOPART_H

// ART includes
namespace fhicl {
  class ParameterSet;
}
#include "cetlib_except/exception.h"

// LArSoft
namespace sim {
  class MCRecoPart;
}

// STL
#include <map>
#include <vector>

namespace sim {

  class MCShowerRecoPart {

  public:
    static const unsigned int kINVALID_UINT;
    static const int kINVALID_INT;

  public:
    /// Default constructor with fhicl parameters
    explicit MCShowerRecoPart(fhicl::ParameterSet const& pset);
    //ClusterMergeAlg(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    /// Main function to read-in data and fill variables in this algorithm to reconstruct MC shower
    void ConstructShower(const MCRecoPart& part_v);

    /**
       Returns a list ot daughter particle index numbers for the specified shower
       with the shower index number as an input
     */
    const std::vector<unsigned int>& ShowerDaughters(const unsigned int shower_id) const
    {
      if (shower_id >= _shower_daughters.size())
        throw cet::exception(__FUNCTION__) << "Invalid shower index!";
      return _shower_daughters.at(shower_id);
    }

    /**
       Returns a list of shower-mother's particle index. Order respects shower index number.
    */
    const std::vector<unsigned int> ShowerMothers() const
    {
      std::vector<unsigned int> mothers(_shower_index.size(), 0);
      for (auto mother_iter = _shower_index.begin(); mother_iter != _shower_index.end();
           ++mother_iter)
        mothers.at((*mother_iter).second) = (*mother_iter).first;
      return mothers;
    }

    //--------------- Particle Information Getters -----------------//

    /**
       Take particle index number and returns shower index number to which this particle belongs.
       Returns -1 if a particle does not belong to any MC shower. Returns kINVALID_INT if input is invalid.
    */
    int ShowerIndex(const unsigned int part_index) const
    {
      if (_shower_id.size() <= part_index) return kINVALID_INT;
      return _shower_id.at(part_index);
    }

  protected:
    /// lots of stdout stream
    bool _debug_mode;

    /// Track index to shower index map
    std::vector<int> _shower_id;

    //
    // shower-indexed-variables
    //
    /// Shower Primary Index ID => Shower Index Map
    std::map<unsigned int, unsigned int> _shower_index;

    /// Shower time-ordered daughters
    std::vector<std::vector<unsigned int>> _shower_daughters;

  }; // class MCShowerRecoPart

} //namespace cluster
#endif
