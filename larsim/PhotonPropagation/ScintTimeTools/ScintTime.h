////////////////////////////////////////////////////////////////////////
// Class:       ScintTime
// Plugin Type: tool
// File:        ScintTime_tool.cc ScintTime.h
// Description:
// Oct. 8 by Mu Wei
////////////////////////////////////////////////////////////////////////
#ifndef ScintTime_H
#define ScintTime_H

namespace CLHEP {
  class HepRandomEngine;
}

namespace phot {
  class ScintTime {
  public:
    virtual ~ScintTime() = default;

    virtual void initRand(CLHEP::HepRandomEngine& engine) = 0;
    virtual void GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine) = 0;
    double GetScintTime() const { return timing; }
    virtual double fastScintTime() = 0;
    virtual double slowScintTime() = 0;

  protected:
    // FIXME: This should be private, with a protected setter.
    //        2022-04-10 CHG
    double timing{0.0};
  };
}

#endif
