
#ifndef TFDH_GSL_FUNCTION_H
#define TFDH_GSL_FUNCTION_H


// simple base class that wraps the gsl_function interface.
// additional parameters to the function are to be held in
// member variables of the derived classes.
class GslFunction {
  public:
  virtual ~GslFunction() = default;
  virtual double f(double x) const = 0;
};


// helper function
// TODO: documentation
double GslFunction_Unpacker(double x, void *params);


#endif // TFDH_GSL_FUNCTION_H
