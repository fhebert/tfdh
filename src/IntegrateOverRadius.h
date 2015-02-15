
#ifndef TFDH_RADIAL_INTEGRAL_H
#define TFDH_RADIAL_INTEGRAL_H

#include <vector>


double integrateOverRadius(const std::vector<double>& r, const std::vector<double>& data);

std::vector<double> accumulateOverRadius(const std::vector<double>& r, const std::vector<double>& data);


#endif // TFDH_RADIAL_INTEGRAL_H

