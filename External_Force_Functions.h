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
    , linear_coeff{linear_coefficient}
    {printf("Force = %fx + %f\n", linear_coeff, constant_coeff);}

    Linear_Force(double linear_coefficient)
    : Linear_Force(linear_coefficient, 0.0) {}

    Linear_Force()
    : Linear_Force(0.0, 0.0) {}

    double force_calculation (double x) override {
        return linear_coeff*x + constant_coeff;
    };
};

class Cubic_Force : public Force_External {
  private:
    double constant_coeff;
    double linear_coeff;
    double quadratic_coeff;
    double cubic_coeff;
  public:
    Cubic_Force(double cubic_coefficient, double quadratic_coefficient
    , double linear_coefficient, double constant_coefficient)
    : constant_coeff{constant_coefficient}
    , linear_coeff{linear_coefficient}
    , quadratic_coeff{quadratic_coefficient}
    , cubic_coeff{cubic_coefficient}
    {printf("Force = %fx^3 + %fx^2 + %fx + %f\n", cubic_coeff, quadratic_coeff, linear_coeff, constant_coeff);}

    double force_calculation (double x) override {
        return cubic_coeff*x*x*x + quadratic_coeff*x*x + linear_coeff*x + constant_coeff;
    };
};

#endif