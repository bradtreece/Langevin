#include <stdio.h>

#ifndef EXTERNAL_FORCE_FUNCTIONS_H
#define EXTERNAL_FORCE_FUNCTIONS_H
class Force_External {
  public:
  Force_External(){}
  virtual double force_calculation (double x)
  {printf("No Function Was Defined For Force Calculation.\n");};
};

class Constant_Force : public Force_External {
  private:
    double constant_value;
  public:
    Constant_Force(double value)
    : constant_value{value}
    {}

    Constant_Force()
    : Constant_Force(0.0) {}

    double force_calculation (double x) override {
        return this->constant_value;
    };
};

class Linear_Force : public Force_External {
  private:
    double constant_coeff;
    double linear_coeff;
  public:
    Linear_Force(double linear_coefficient, double constant_coefficient)
    : constant_coeff{constant_coefficient}
    , linear_coeff{linear_coeff}
    {}

    Linear_Force(double linear_coefficient)
    : Linear_Force(linear_coefficient, 0.0) {}

    Linear_Force()
    : Linear_Force(0.0, 0.0) {}

    double force_calculation (double x) override {
        return linear_coeff*x + constant_coeff;
    };
};

#endif