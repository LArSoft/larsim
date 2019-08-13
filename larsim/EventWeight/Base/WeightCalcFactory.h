#ifndef _WEIGHTCALCFACTORY_H_
#define _WEIGHTCALCFACTORY_H_

#include <map>
#include <string>

namespace evwgh {
  class WeightCalc;
  class WeightCalcCreator;

  class WeightCalcFactory
  {
  public:
    static WeightCalc* Create(const std::string& classname);
    static void Register(const std::string& wghcalcname,
			 WeightCalcCreator* creator);

  private:
    static std::map<std::string, WeightCalcCreator*>& GetTable();
  };
}

#endif // _WEIGHTCALCFACTORY_H_
