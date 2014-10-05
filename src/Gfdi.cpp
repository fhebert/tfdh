
#include "Gfdi.h"

#include <array>
#include <cassert>
#include <cmath>


// helper functions and data
namespace {

  // numerous constants used in the analytic approximations
  //const std::array<double, 5> x = {{
  const auto x = std::array<double, 5> {{
    7.265351e-2, 0.2694608, 0.533122, 0.7868801, 0.9569313}};

  const auto xi = std::array<double, 5> {{
    0.26356032, 1.4134031, 3.5964258, 7.0858100, 12.640801}};

  const auto h = std::array<double, 5> {{
    3.818735e-2, 0.1256732, 0.1986308, 0.1976334, 0.1065420}};

  const auto v = std::array<double, 5> {{
    0.29505869, 0.32064856, 7.3915570e-2, 3.6087389e-3, 2.3369894e-5}};

  const auto c = std::array<std::array<double, 5>, 3> {{
    {{0.37045057, 0.41258437, 9.777982e-2, 5.3734153e-3, 3.8746281e-5}},
    {{0.39603109, 0.69468795, 0.22322760, 1.5262934e-2, 1.3081939e-4}},
    {{0.76934619, 1.7891437, 0.70754974, 5.6755672e-2, 5.5571480e-4}}}};

  const auto khi = std::array<std::array<double, 5>, 3> {{
    {{0.43139881, 1.7597537, 4.1044654, 7.7467038, 13.457678}},
    {{0.81763176, 2.4723339, 5.1160061, 9.0441465, 15.049882}},
    {{1.2558461, 3.2070406, 6.1239082, 10.316126, 16.597079}}}};


  inline double cube(const double x) {
    return x*x*x;
  }

  // helper function for GFDI work
  double gfdi_helper(const int k, const double chi, const double tau,
      const double r) {
    if (chi*tau < 1.e-4 && chi > 0.0) {
      return pow(chi, k+3./2)/(k+3./2);
    }
    else if (k==0) {
      return (chi + 1/tau)*r/2
        - pow(2*tau, -3./2) * log(1 + tau*chi + sqrt(2*tau)*r);
    }
    else if (k==1) {
      return (2./3*cube(r) - gfdi_helper(0, chi, tau, r)) / tau;
    }
    else if (k==2) {
      return (2*chi*cube(r) - 5*gfdi_helper(1, chi, tau, r)) / (4*tau);
    }
    else {
      return 0.0;
    }
  }

} // end anonymous namespace



double gfdi(const GFDI order, const double chi, const double tau) {
  // TODO: something more like a SpEC require?
  assert(tau <= 100. && "GFDI: outside of known convergence region");

  const int k = static_cast<int>(order);

  if (chi <= 0.6) {
    double value = 0;
    for (int i=1; i<=5; ++i) {
      value += c[k][i-1] * sqrt(1 + khi[k][i-1]*tau/2) /
        (exp(-khi[k][i-1]) + exp(-chi));
    }
    return value;
  }
  else if (chi < 14.0) {
    double value = 0;
    for (int i=1; i<=5; ++i) {
      value += h[i-1] * pow(x[i-1], k) * pow(chi, k+3./2)
        * sqrt(1 + chi*x[i-1]*tau/2) / (1 + exp(chi*(x[i-1] - 1)))
        + v[i-1] * pow(xi[i-1] + chi, k+1./2) * sqrt(1 + (xi[i-1] + chi)*tau/2);
    }
    return value;
  }
  else {
    const double r = sqrt(chi*(1 + chi*tau/2));
    return gfdi_helper(k, chi, tau, r)
      + M_PI*M_PI/6. * pow(chi, k) * (k + 1./2 + (k+1)*chi*tau/2) / r;
  }
}
