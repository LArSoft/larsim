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

// Random number engine
#include "CLHEP/Random/RandFlat.h"

#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"

#include <memory>

namespace fhicl {
  class ParameterSet;
}

namespace phot {
  class ScintTimeLAr : public ScintTime {
  public:
    explicit ScintTimeLAr(fhicl::ParameterSet const& pset);
    void initRand(CLHEP::HepRandomEngine& engine);
    void GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine);
    double fastScintTime();
    double slowScintTime();

  private:
    int LogLevel;
    std::unique_ptr<CLHEP::RandFlat> fUniformGen;
    // parameters for the shape of argon scintillation light time distribution
    const double SRTime; // PureLAr: rising time of slow LAr scintillation;
    const double SDTime; // PureLAr: decay time of slow LAr scintillation;
    const double FRTime; // PureLAr: rising time of fast LAr scintillation;
    const double FDTime; // PureLAr: decay time of fast LAr scintillation;
    const bool fNoFastRisingTime, fNoSlowRisingTime;

    // general functions
    double single_exp(double t, double tau2) const;
    double bi_exp(double t, double tau1, double tau2) const;
    double with_rising_time(double tau1, double tau2);
  };
}
#endif
