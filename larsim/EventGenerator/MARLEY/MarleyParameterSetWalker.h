//////////////////////////////////////////////////////////////////////////////
/// \file MarleyParameterSetWalker.h
/// \brief Concrete fhicl::ParameterSetWalker that converts a
/// fhicl::ParameterSet into a marley::JSON object
///
/// \author Steven Gardiner <gardiner@fnal.gov>
//////////////////////////////////////////////////////////////////////////////

#ifndef LARSIM_MARLEY_PARAMETERSET_WALKER_H
#define LARSIM_MARLEY_PARAMETERSET_WALKER_H

// standard library includes
#include <sstream>
#include <string>
#include <vector>

// framework includes
#include "fhiclcpp/ParameterSetWalker.h"
#include "fhiclcpp/detail/printing_helpers.h"

// MARLEY includes
#include "marley/JSON.hh"

namespace evgen {

  class MarleyParameterSetWalker : public fhicl::ParameterSetWalker {

  public:
    inline const marley::JSON& get_json() const { return full_json_; }
    inline marley::JSON& get_json() { return full_json_; }

  private:
    inline void enter_table(const key_t& key, const any_t&) override
    {
      auto& json = json_refs_.back().get();
      json_refs_.emplace_back(json[key] = marley::JSON::make(marley::JSON::DataType::Object));
    }

    inline void enter_sequence(const key_t& key, const any_t&) override
    {
      auto& json = json_refs_.back().get();
      in_seq_ = true;
      seq_index_ = 0u;
      json_refs_.emplace_back(json[key] = marley::JSON::make(marley::JSON::DataType::Array));
    }

    inline void atom(const key_t& key, const any_t& any) override
    {
      auto& json = json_refs_.back().get();
      // Convert the atom to a string. Add an extra space to keep some
      // of MARLEY's JSON parsing routines (which rely on looking ahead
      // in some cases) happy
      std::string atom_val = fhicl::detail::atom::value(any) + ' ';
      marley::JSON json_atom;
      // Hard-coded value taken from fhiclcpp/detail/printing_helpers.cc
      if (atom_val != "@nil") {
        std::istringstream iss(atom_val);
        // Utility function defined in marley/JSON.hh
        json_atom = marley::parse_next(iss);
      }
      if (in_seq_) { json[seq_index_++] = json_atom; }
      else {
        json[key] = json_atom;
      }
    }

    inline void exit_table(const key_t&, const any_t&) override { json_refs_.pop_back(); }

    inline void exit_sequence(const key_t&, const any_t&) override
    {
      json_refs_.pop_back();
      in_seq_ = false;
    }

    /// Owned JSON object used to store the converted FHiCL parameters
    marley::JSON full_json_;

    /// References to the owned JSON object or a sub-object thereof
    std::vector<std::reference_wrapper<marley::JSON>> json_refs_ = {full_json_};

    unsigned seq_index_ = 0u;
    bool in_seq_ = false;
  };

}

#endif // LARSIM_MARLEY_PARAMETERSET_WALKER_H
