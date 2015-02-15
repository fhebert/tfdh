
#ifndef TFDH_TFDH_SOLUTION_H
#define TFDH_TFDH_SOLUTION_H

#include <cassert>
#include <vector>


struct TfdhSolution {
  const std::vector<double> r;
  const std::vector<double> phi;

  TfdhSolution(const std::vector<double>& rs, const std::vector<double>& phis)
    : r(rs), phi(phis)
  {
    assert(r.size()==phi.size());
  }
};


#endif // TFDH_TFDH_SOLUTION_H
