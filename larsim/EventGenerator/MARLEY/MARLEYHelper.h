//////////////////////////////////////////////////////////////////////////////
/// \file MARLEYHelper.h
/// \brief LArSoft interface to the MARLEY (Model of Argon Reaction Low Energy
/// Yields) supernova neutrino event generator
///
/// \author Steven Gardiner <sjgardiner@ucdavis.edu>
//////////////////////////////////////////////////////////////////////////////

#ifndef LARSIM_ALGORITHMS_MARLEYGENERATOR_H
#define LARSIM_ALGORITHMS_MARLEYGENERATOR_H

// standard library includes
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// framework includes
namespace fhicl {
  class ParameterSet;
}

// art extensions
namespace rndm {
  class NuRandomService;
}

// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TLorentzVector.h"

// MARLEY includes
#include "marley/Generator.hh"

namespace marley {
  class Event;
  class JSON;
  class Particle;
}

namespace evgen {

  class MARLEYHelper {

  public:
    MARLEYHelper(const fhicl::ParameterSet& pset,
                 rndm::NuRandomService& rand_service,
                 const std::string& generator_name);

    void reconfigure(const fhicl::ParameterSet& pset);

    // If a non-null marley::Event* is supplied, the marley::Event
    // object corresponding to the generated MCTruth object is loaded
    // into the target of the pointer.
    simb::MCTruth create_MCTruth(const TLorentzVector& vtx_pos,
                                 marley::Event* marley_event = nullptr);

    marley::Generator& get_generator() { return *fMarleyGenerator; }
    const marley::Generator& get_generator() const { return *fMarleyGenerator; }

    std::string find_file(const std::string& fileName, const std::string& fileType);

  protected:
    void add_marley_particles(simb::MCTruth& truth,
                              const std::vector<marley::Particle*>& particles,
                              const TLorentzVector& vtx_pos,
                              bool track);

    void load_full_paths_into_json(marley::JSON& json,
                                   const std::string& array_name,
                                   bool missing_ok = false);

    std::unique_ptr<marley::Generator> fMarleyGenerator;

    // name to use for this instance of MARLEYHelper
    std::string fHelperName;

    // string stream used to capture logger output from MARLEY
    // and redirect it to the LArSoft logger
    std::stringstream fMarleyLogStream;

    // Loads ROOT dictionaries for the MARLEY Event and Particle classes.
    // This allows a module to write the generated events to a TTree.
    void load_marley_dictionaries();

  }; // class evgen::MARLEYHelper

} // namespace evgen

#endif // LARSIM_ALGORITHMS_MARLEYGENERATOR_H
