#include <stdio.h>

class MODEL {
public:
  MODEL(){};
  virtual void model_call (double *t, double *x, int length, double *params)
  {printf("No Call To This Model Was Defined.\n");};
};

class Linear_Model : public MODEL {
public:
  void model_call(double *t, double *x, int length, double *params);
};