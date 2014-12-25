
#include "TfdhFunctions.h"

#include "GslWrappers.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "PhysicalConstants.h"
#include "RadialFunction.h"

#include <cmath>

#include <iostream>

namespace {

  class FermiDiracDistribution : public GslFunction {
    const double xi;
    const PlasmaState p;
    public:
    FermiDiracDistribution(const double xi, const PlasmaState& p)
      : xi(xi), p(p) {}
    double f(const double x) const {
      const double val = (1 + p.tau*x) * sqrt(x + p.tau * x*x/2.0) / (1.0 + exp(x-p.chi-xi));
      //std::cout << " FermiDirac(xi=" << xi << ") computed at x=" << x << ", got val=" << val << std::endl;
      //std::cout << "   x-chi-xi = " << x-p.chi-xi << ",  exp(x-chi-xi) = " << exp(x-p.chi-xi) << std::endl;
      return val;
    }
  };

  double boundFermiDiracIntegral(const double xi, const PlasmaState& p) {
    FermiDiracDistribution fd(xi, p);
    const double eps = 1.e-6;
    return gslQuadratureAG(fd, 0, xi, eps, eps);
  }

}

double TFDH::neBound(const double phi, const PlasmaState& p)
{
  const double xi = phi/p.kt;
  // TODO: check these conditions!
  if (xi < 0) {
    return 0;
  } else if (xi > 1.e4*(1+fabs(p.chi))) {
    // xi very large => definite integral approaches indefinite integral
    return Plasma::ne(phi, p);
  } else {
    const double& NePrefactor = PhysicalConstantsCGS::NePrefactor;
    return NePrefactor * pow(p.kt, 1.5) * boundFermiDiracIntegral(xi, p);
  }
}

RadialFunction TFDH::neBound(const RadialFunction& tfdh, const PlasmaState& p)
{
  /*
  const int place = 1;
  const double xi = tfdh.data[place]/tfdh.radii[place]/p.kt;
  const FermiDiracDistribution fd(xi,p);
  const int len = 1000;
  std::vector<double> xs(len), fs(len);
  for (size_t i=0; i<len; ++i) {
    const double x = i/(len-1.0) * xi;
    xs[i] = x;
    fs[i] = fd.f(x);
  }
  const RadialFunction test(xs, fs);
  writeToFile(test, "test.dat");
  return test;
  */
  std::vector<double> nebs(tfdh.data.size());
  for (size_t i=0; i<tfdh.radii.size(); ++i) {
    //std::cout << "computing neB..." << std::endl;
    nebs[i] = neBound(tfdh.data[i], p);
    //std::cout << "  success" << std::endl;
  }
  return RadialFunction(tfdh.radii, nebs);
}


