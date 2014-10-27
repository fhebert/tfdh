
#ifndef TFDH_INTEGRATE_TFDH_HPP
#define TFDH_INTEGRATE_TFDH_HPP

struct Element;
struct PlasmaState;
struct RadialFunction;

namespace TFDH {
  RadialFunction solve(const Element& e, const PlasmaState& p);
  RadialFunction integrateODE(const Element& e, const PlasmaState& p,
      double r_init, double r_final, double dv0);
  double rootfindPotential(const Element& e, const PlasmaState& p,
      double r_init, double r_final);
};

#endif // TFDH_INTEGRATE_TFDH_HPP
