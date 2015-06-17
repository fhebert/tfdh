
#ifndef TFDH_TFDH_SOLUTION_H
#define TFDH_TFDH_SOLUTION_H

#include "GslWrappers.h"

#include <cassert>
#include <vector>


class TfdhSolution : public GSL::FunctionObject {
  public:
    TfdhSolution(const std::vector<double>& rs, const std::vector<double>& phis)
      : r(rs), phi(phis), spline(r, phi)
    {
      assert(r.size()==phi.size());
      for (size_t i=0; i<r.size()-1; ++i)
        assert(r[i] < r[i+1] and "need monotonically increasing r");
    }

    double operator()(double radius) const override {return spline.eval(radius);}

  public:
    const std::vector<double> r;
    const std::vector<double> phi;

  private:
    GSL::Spline spline;
};


#endif // TFDH_TFDH_SOLUTION_H
