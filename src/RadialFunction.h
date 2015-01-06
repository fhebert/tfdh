
#ifndef TFDH_RADIAL_FUNCTION_H
#define TFDH_RADIAL_FUNCTION_H

#include <cassert>
#include <string>
#include <vector>


struct RadialFunction {
  const std::vector<double> radii;
  const std::vector<double> data;

  RadialFunction(const std::vector<double>& rs, const std::vector<double>& fs)
    : radii(rs), data(fs)
  {
    assert(rs.size()==fs.size());
  }
};


int writeToFile(const RadialFunction& rf, const std::string& filename);


#endif // TFDH_RADIAL_FUNCTION_H

