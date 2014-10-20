
#ifndef TFDH_INTEGRATE_TFDH_HPP
#define TFDH_INTEGRATE_TFDH_HPP

struct Element;
struct PlasmaState;
struct RadialFunction;

namespace TFDH {
  RadialFunction integrate(const Element& e, const PlasmaState& p);
};

#endif // TFDH_INTEGRATE_TFDH_HPP
