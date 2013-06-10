#include <cmath>
#include "initialisers.h"

const double pi = 3.14159;

double sin_function(double x, double y)
{
  double value = (2*pi*pi*(sin(x*pi) + sin(y*pi)));
  return value;
}

double zero(double x, double y)
{
  return 0.0;
}

double sin_function_2(double x, double y)
{
  double value = ((sin(x*pi) + sin(y*pi)));
  return value;
}

double sinh_function(double x, double y)
{
  double value = ((sin(x*pi)*sinh(y*pi)));
  return value;
}

double constant(double x, double y)
{
  double value = 1.0;
  return value;
}

double constant2(double x, double y)
{
  return x + y;
}

double ramp(double x, double y)
{
  return x;
}

double omega(double x, double y)
{
  double delta = 0.01;
  if((std::abs(0.5 - x) < delta) && (std::abs(0.5 - y) < delta))
  {
    return 10e6;
  }
  else
  {
    return 1;
  }
}

double pi_function(double x, double y)
{
  return pi*pi;
}