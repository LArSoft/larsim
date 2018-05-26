#include <map>
#include <string>
#include "canvas/Persistency/Common/Wrapper.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

template class art::Wrapper<evwgh::MCEventWeight>;
template class art::Wrapper<std::vector<evwgh::MCEventWeight> >;
template class std::map<std::string, double>;
template class std::map<std::string, std::vector<double> >;
