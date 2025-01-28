#ifndef _WEIGHTCALCCREATOR_H_
#define _WEIGHTCALCCREATOR_H_

#include <iostream>
#include <string>

namespace evwgh {
  class WeightCalc;

  class WeightCalcCreator {
  public:
    WeightCalcCreator(const std::string& classname);
    virtual ~WeightCalcCreator() = default;

    virtual WeightCalc* Create() = 0;
  };

  template <class T>
  class WeightCalcImpl : public WeightCalcCreator {
  public:
    WeightCalcImpl(const std::string& classname) : WeightCalcCreator(classname)
    {
      std::cout << __PRETTY_FUNCTION__ << " called." << std::endl;
    }
    virtual ~WeightCalcImpl() {}
    virtual WeightCalc* Create() { return new T; }
  };

#define DECLARE_WEIGHTCALC(wghcalc)                                                                \
private:                                                                                           \
  static const WeightCalcImpl<wghcalc> creator;

#define REGISTER_WEIGHTCALC(wghcalc) const WeightCalcImpl<wghcalc> wghcalc::creator(#wghcalc);
}

#endif // _WEIGHTCALCFACTORY_H_
