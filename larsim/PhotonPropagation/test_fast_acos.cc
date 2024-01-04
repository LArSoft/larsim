#include "larsim/PhotonPropagation/PhotonPropagationUtils.h"

#include <cmath>
#include <iomanip>
#include <iostream>


using phot::fast_acos;

int main()
{
  std::cout << "fast\tstd\n";
  std::cout << std::setprecision(17);
  for (int i = 0; i <= 100; ++i)
  {
    double x = static_cast<double>(i)/100.0;
    double f_x = fast_acos(x);
    double y_x = std::acos(x);
    std::cout << f_x << '\t' << y_x << '\n';
  }
}
