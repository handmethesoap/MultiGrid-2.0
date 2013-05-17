#include "grid.h"
#include <iostream>
  
Grid:: Grid(int levels)
{
  int length_v = 0;
  int length_f = 0; 
  int temp = 0;
  
  for(int i = 1; i <= levels; ++i)
  {
    temp = (1 << i);
    length_v += (temp + 1)*(temp + 1);
    length_f += (temp - 1)*(temp - 1);
  }
  
  _v = new double[length_v];
  _f = new double[length_f];
  
  std::cout << length_v << ", " << length_f << std::endl;
}

Grid:: ~Grid()
{
  delete[] _v;
  delete[] _f;
}