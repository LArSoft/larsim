#ifndef _MCEVENTWEIGHT_H_
#define _MCEVENTWEIGHT_H_

#include <vector>
#include <string>

namespace evwgh {
  struct MCEventWeight
  {
    std::map<std::string, std::vector<double> > fWeight;
  };
}
#endif //_MCEVENTWEIGHT_H_
