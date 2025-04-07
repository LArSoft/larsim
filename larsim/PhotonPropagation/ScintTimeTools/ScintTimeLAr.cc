////////////////////////////////////////////////////////////////////////
// Class:       ScintTimeLAr
// Plugin Type: tool
// File:        ScintTimeLAr_tool.cc ScintTimeLAr.h
// Description:
// Generate a random number for timing of LAr scintillation
// Oct. 8 by Mu Wei
////////////////////////////////////////////////////////////////////////
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTimeLAr.h"

#include "fhiclcpp/ParameterSet.h"

namespace phot {
  //......................................................................
  ScintTimeLAr::ScintTimeLAr(fhicl::ParameterSet const& pset)
    : LogLevel{pset.get<int>("LogLevel")}
    , SRTime{pset.get<double>("SlowRisingTime", 0.)}
    , SDTime{pset.get<double>("SlowDecayTime", 0.)}
    , FRTime{pset.get<double>("FastRisingTime", 0.)}
    , FDTime{pset.get<double>("FastDecayTime", 0.)}
    , fNoFastRisingTime(pset.get<bool>("NoFastRisingTime", false))
    , fNoSlowRisingTime(pset.get<bool>("NoSlowRisingTime", false))
  {
    if (LogLevel >= 1) {
      std::cout << "ScintTimeLAr Tool configure:" << std::endl;
      std::cout << "Fast rising time: " << FRTime << ", Fast decay time: " << FDTime
                << ", Slow rising time: " << SRTime << ", Slow decay time: " << SDTime << std::endl;
    }
  }

  void ScintTimeLAr::initRand(CLHEP::HepRandomEngine& engine)
  {
    fUniformGen = std::make_unique<CLHEP::RandFlat>(engine);
  }

  //......................................................................
  double ScintTimeLAr::single_exp(double t, double tau2) const
  {
    return std::exp((-1.0 * t) / tau2) / tau2;
  }

  //......................................................................
  double ScintTimeLAr::bi_exp(double t, double tau1, double tau2) const
  {
    return (((std::exp((-1.0 * t) / tau2) * (1.0 - std::exp((-1.0 * t) / tau1))) / tau2) / tau2) *
           (tau1 + tau2);
  }

  //......................................................................
  double ScintTimeLAr::with_rising_time(double tau1, double tau2)
  {
    // TODO: This method is very inefficient.
    // If we care about simulating the rising time, it would be better
    // to simulate random numbers according a general distribution
    double d = (tau1 + tau2) / tau2;
    int ncalls = 0;
    std::cout << "Calling with_rising_time" << std::endl;
    while (1) {
      ++ncalls;
      double ran1 = fUniformGen->fire();
      double ran2 = fUniformGen->fire();
      double t = -tau2 * std::log(1 - ran1);
      double g = d * single_exp(t, tau2);
      if (ran2 <= bi_exp(t, tau1, tau2) / g) {
        std::cout << "sim'd rise time after " << ncalls << " calls" << std::endl;
        return t;
      }
    }
  }

  //......................................................................
  // Returns the time within the time distribution of the
  // scintillation process, when the photon was created.
  // Scintillation light has an exponential decay which is given by
  // the decay time, tau2, and an exponential increase, which here is
  // given by the rise time, tau1.
  void ScintTimeLAr::GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine)
  {
    double tau1;
    double tau2;

    if (is_fast) {
      tau1 = FRTime;
      tau2 = FDTime;
    }
    else {
      tau1 = SRTime;
      tau2 = SDTime;
    }

    CLHEP::RandFlat randflatscinttime{engine};

    if ((tau1 < 1e-8) || (tau1 == -1.0)) {
      timing = -tau2 * std::log(randflatscinttime());
      return;
    }

    //ran1, ran2 = random numbers for the algorithm
    while (1) {
      auto ran1 = randflatscinttime();
      auto ran2 = randflatscinttime();
      auto d = (tau1 + tau2) / tau2;
      auto t = -tau2 * std::log(1 - ran1);
      auto g = d * single_exp(t, tau2);
      if (ran2 <= bi_exp(t, tau1, tau2) / g) {
        timing = t;
        return;
      }
    }
  }

  inline double ScintTimeLAr::fastScintTime()
  {
    if (fNoFastRisingTime) { return -FDTime * std::log(fUniformGen->fire()); }
    return with_rising_time(FRTime, FDTime);
  }

  inline double ScintTimeLAr::slowScintTime()
  {
    if (fNoSlowRisingTime) {
      return -SDTime * std::log(fUniformGen->fire());
    }
    return with_rising_time(SRTime, SDTime);
  }

}
