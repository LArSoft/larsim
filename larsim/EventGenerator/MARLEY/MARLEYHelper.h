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
#include "TTree.h"

// HepMC3 includes
#include "HepMC3/Data/GenEventData.h"

// MARLEY includes
#include "marley/Generator.hh"

namespace marley {
  class JSON;
}

namespace evgen {

  class MARLEYHelper {

  public:
    // If a non-nullptr TTree* is supplied, then the run information
    // (HepMC3::GenRunInfoData) and generated events (HepMC3::GenEventData)
    // will be dumped in a format similar to standalone "marley generate"
    MARLEYHelper(const fhicl::ParameterSet& pset,
                 rndm::NuRandomService& rand_service,
                 const std::string& generator_name,
                 TTree* dump_tree);

    void reconfigure(const fhicl::ParameterSet& pset);

    // object corresponding to the generated MCTruth object is loaded
    // into the target of the pointer.
    simb::MCTruth create_MCTruth(const TLorentzVector& vtx_pos);

    marley::Generator& get_generator() { return *fMarleyGenerator; }
    const marley::Generator& get_generator() const { return *fMarleyGenerator; }

    std::string find_file(const std::string& fileName, const std::string& fileType);

  protected:
    void add_marley_particles(simb::MCTruth& truth,
                              const std::vector<std::shared_ptr<HepMC3::GenParticle>>& particles,
                              const TLorentzVector& vtx_pos,
                              double conv_factor,
                              bool track);

    void clear_event_data();

    void load_full_paths_into_json(marley::JSON& json,
                                   const std::string& array_name,
                                   bool missing_ok = false);

    std::unique_ptr<marley::Generator> fMarleyGenerator;

    // name to use for this instance of MARLEYHelper
    std::string fHelperName;

    // string stream used to capture logger output from MARLEY
    // and redirect it to the LArSoft logger
    std::stringstream fMarleyLogStream;

    // TTree to use for optionally dumping event and run information
    TTree* fDumpTree = nullptr;

    // Storage used to help with dumping
    std::shared_ptr<HepMC3::GenRunInfo> fRunInfo;
    std::unique_ptr<HepMC3::GenEventData> fEventData;
    HepMC3::GenEventData* fEventDataPtr = nullptr;

  }; // class evgen::MARLEYHelper

} // namespace evgen

#endif // LARSIM_ALGORITHMS_MARLEYGENERATOR_H
