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
#include <array>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// framework includes
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/detail/ParameterArgumentTypes.h"
#include "fhiclcpp/types/detail/ParameterBase.h"
#include "fhiclcpp/types/detail/TableBase.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT includes
#include "TLorentzVector.h"

// MARLEY includes
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/Particle.hh"

namespace evgen {

  class MARLEYHelper {

    public:

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      /// Collection of configuration parameters that will be
      /// forwarded to MARLEY and used to define the neutrino
      /// source
      struct Source_Config {
        fhicl::Atom<std::string> type_ {
          Name("type"),
          Comment("Type of neutrino source for MARLEY to use")
        };

        fhicl::Atom<std::string> neutrino_ {
          Name("neutrino"),
          Comment("Kind of neutrino (flavor, matter/antimatter) that the"
            " neutrino source produces")
        };

        fhicl::Atom<double> Emin_ {
          Name("Emin"),
          Comment("Minimum energy (MeV) of the neutrinos produced by this"
            " source"),
          [this]() -> bool {
            auto type = type_();
            return (type == "fermi-dirac") || (type == "beta-fit" );
          }
        };

        fhicl::Atom<double> Emax_ {
          Name("Emax"),
          Comment("Maximum energy (MeV) of the neutrinos produced by this"
            " source"),
          [this]() -> bool {
            auto type = type_();
            return (type == "fermi-dirac") || (type == "beta-fit" )
              || (type == "histogram");
          }
        };

        fhicl::Atom<double> temperature_ {
          Name("temperature"),
          Comment("Effective temperature for the Fermi-Dirac distribution"),
          [this]() -> bool {
            auto type = type_();
            return (type == "fermi-dirac");
          }
        };

        fhicl::OptionalAtom<double> eta_ {
          Name("eta"),
          Comment("Pinching parameter for the Fermi-Dirac distribution"),
          [this]() -> bool {
            auto type = type_();
            return (type == "fermi-dirac");
          }
        };

        fhicl::Atom<double> energy_ {
          Name("energy"),
          Comment("Energy (MeV) of the neutrinos produced by a monoenergetic"
            " source"),
          [this]() -> bool {
            auto type = type_();
            return (type == "monoenergetic");
          }
        };

        fhicl::Atom<double> Emean_ {
          Name("Emean"),
          Comment("Mean energy (MeV) of the neutrinos produced by a beta-fit"
            " source"),
          [this]() -> bool {
            auto type = type_();
            return (type == "beta-fit");
          }
        };

        fhicl::OptionalAtom<double> beta_ {
          Name("beta"),
          Comment("Pinching parameter for a beta-fit source"),
          [this]() -> bool {
            auto type = type_();
            return (type == "beta-fit");
          }
        };

        fhicl::Sequence<double> E_bin_lefts_ {
          Name("E_bin_lefts"),
          Comment("Left edges for each energy bin in the histogram"),
          [this]() -> bool {
            auto type = type_();
            return type == "histogram";
          }
        };

        fhicl::Sequence<double> weights_ {
          Name("weights"),
          Comment("Weights for each energy bin in the histogram"),
          [this]() -> bool {
            auto type = type_();
            return type == "histogram";
          }
        };

        fhicl::Sequence<double> energies_ {
          Name("energies"),
          Comment("Energies (MeV) for each grid point"),
          [this]() -> bool {
            auto type = type_();
            return type == "grid";
          }
        };

        fhicl::Sequence<double> prob_densities_ {
          Name("prob_densities"),
          Comment("Probability densities for each grid point"),
          [this]() -> bool {
            auto type = type_();
            return type == "grid";
          }
        };

        fhicl::Atom<std::string> rule_ {
          Name("rule"),
          Comment("Interpolation rule for computing probability densities"
            " between grid points. Allowed values include \"linlin\","
            " \"linlog\", \"loglin\", \"loglog\", and \"constant\""),
          [this]() -> bool {
            auto type = type_();
            return type == "grid";
          }
        };

        fhicl::Atom<std::string> tfile_ {
          Name("tfile"),
          Comment("Name of the ROOT file that contains a TH1 or TGraph neutrino"
            " source spectrum"),
          [this]() -> bool {
            auto type = type_();
            return (type == "th1") || (type == "tgraph");
          }
        };

        fhicl::Atom<std::string> namecycle_ {
          Name("namecycle"),
          Comment("Namecycle of the ROOT TH1 or TGraph to use"),
          [this]() -> bool {
            auto type = type_();
            return (type == "th1") || (type == "tgraph");
          }
        };

        fhicl::Atom<bool> weight_flux_ {
          Name("weight_flux"),
          Comment("Whether to weight the input flux by the reaction cross section(s)"),
          true
        };

      }; // struct Source_Config

      /// Collection of configuration parameters that will be
      /// forwarded to MARLEY
      struct Config {

        fhicl::OptionalAtom<std::string> seed_ {
          Name("seed"),
          Comment("Seed used for the MARLEY generator")
        };

        fhicl::Sequence<double, 3> direction_ {
          Name("direction"),
          Comment("3-vector that points in the direction of motion of the"
            " incident neutrinos"),
	  // for c2: need double braces
          std::array<double, 3> {{ 0., 0., 1. }} // default value
        };

        fhicl::Sequence<std::string> reactions_ {
          Name("reactions"),
          Comment("List of matrix element data files to use to define reactions"
            " in MARLEY")
        };

        fhicl::Sequence<std::string> structure_ {
          Name("structure"),
          Comment("List of TALYS format nuclear structure data files to use"
            " in MARLEY")
        };

        fhicl::Table<Source_Config> source_ {
          Name("source"),
          Comment("Neutrino source configuration")
        };

      }; // struct Config

      // Configuration-checking constructors
      MARLEYHelper(const fhicl::Table<Config>& conf,
        rndm::NuRandomService& rand_service,
        const std::string& generator_name);

      MARLEYHelper(const fhicl::ParameterSet& pset,
        rndm::NuRandomService& rand_service, const std::string& generator_name)
        : MARLEYHelper(fhicl::Table<Config>(pset, {}),
        rand_service, generator_name) {}

      void reconfigure(const fhicl::Table<Config>& conf);

      // If a non-null marley::Event* is supplied, the marley::Event
      // object corresponding to the generated MCTruth object is loaded
      // into the target of the pointer.
      simb::MCTruth create_MCTruth(const TLorentzVector& vtx_pos,
        marley::Event* marley_event = nullptr);

      marley::Generator& get_generator() { return *fMarleyGenerator; }
      const marley::Generator& get_generator() const
        { return *fMarleyGenerator; }

      std::string find_file(const std::string& fileName,
        const std::string& fileType);

    protected:

      template <typename T> marley::JSON fhicl_atom_to_json(
        const fhicl::Atom<T>* atom)
      {
        return marley::JSON(atom->operator()());
      }

      template <typename T> marley::JSON fhicl_optional_atom_to_json(
        const fhicl::OptionalAtom<T>* opt_atom)
      {
        T value;
        if (opt_atom->operator()(value)) return marley::JSON(value);
        // Return a null JSON object if the fhicl::OptionalAtom wasn't used
        else return marley::JSON();
      }

      template <typename S> marley::JSON fhicl_sequence_to_json(
        const S* sequence)
      {
        marley::JSON array = marley::JSON::make(
          marley::JSON::DataType::Array);

        for (size_t i = 0; i < sequence->size(); ++i) {
          array.append(sequence->operator()(i));
        }
        return array;
      }

      marley::JSON fhicl_parameter_to_json(
        const fhicl::detail::ParameterBase* par);

      void add_marley_particles(simb::MCTruth& truth,
        const std::vector<marley::Particle*>& particles,
        const TLorentzVector& vtx_pos, bool track);

      void load_full_paths_into_json(marley::JSON& json,
        const std::string& array_name);

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
