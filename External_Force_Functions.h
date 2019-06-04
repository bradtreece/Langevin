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

#endif