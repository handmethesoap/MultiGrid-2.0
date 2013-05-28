#include "stdafx.h"
#include <cmath>
#include "initialisers.h"

const double pi = 3.14159;

double sin_function(double x, double y)
{
  double value = (2*pi*pi*(sin(x*pi) + sin(y*pi)));
  return value;
}

double constant(double x, double y)
{
  double value = 1.0;
  return value;
}

double constant2(double x, double y)
{
	static int number = 32;
	int temp;
	temp = number;
	number = 0;
	return temp;
}