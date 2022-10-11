////////////////////////////////////////////////////////////////////////
// Class:       ScintTimeLAr
// Plugin Type: tool
// File:        ScintTimeLAr_tool.cc ScintTimeLAr.h
// Description:
// Generate a random number for timing of LAr scintillation
// Oct. 8 by Mu Wei
////////////////////////////////////////////////////////////////////////
#ifndef ScintTimeLAr_H
#define ScintTimeLAr_H

#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"

namespace fhicl {
  class ParameterSet;
}

namespace phot {
  class ScintTimeLAr : public ScintTime {
  public:
    explicit ScintTimeLAr(fhicl::ParameterSet const& pset);
    void GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine);

  private:
    int LogLevel;

    // parameters for the shape of argon scinitllation light time distribution
    double SRTime; // PureLAr: rising time of slow LAr scinitllation;
    double SDTime; // PureLAr: decay time of slow LAr scintillation;
    double FRTime; // PureLAr: rising time of fast LAr scinitllation;
    double FDTime; // PureLAr: decay time of fast LAr scintillation;

    // general functions
    double single_exp(double t, double tau2) const;
    double bi_exp(double t, double tau1, double tau2) const;
  };
}
#endif
