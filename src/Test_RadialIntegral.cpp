
#include "RadialFunction.h"
#include "RadialIntegral.h"

#include <cmath>
#include <string>
#include <iostream>
#include <vector>

int main() {
  const int numpts = 200;
  std::vector<double> r(numpts), f(numpts);
  for (int i=0; i<numpts; ++i) {
    const double ri = i/(numpts-1.0);
    r[i] = ri;
    f[i] = ri*ri - ri;
  }
  const RadialFunction fr(r, f);
  const double integral = RadialIntegral(fr);
  std::cout << "integral value = " << integral << ", vs expected = " << -M_PI/5.0 << "\n";

  return 0;
}
