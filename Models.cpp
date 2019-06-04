#include "Models.h"

void Constant_Model::model_call(double *t, double *x, int length, double *params)
{
  for (int i=0; i<length; i++)
  {
    x[i] = params[0];
  }
}

void Linear_Model::model_call(double *t, double *x, int length, double *params)
{
  for (int i=0; i<length; i++)
  {
      x[i] = params[1]*t[i] + params[0];
  }
}