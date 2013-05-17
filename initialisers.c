#include <cmath>
#include "initialisers.h"

const double pi = 3.14159;

double sin_function(double x, double y)
{
  double value = (2*pi*pi*(sin(x*pi) + sin(y*pi)));
  return value;
}
