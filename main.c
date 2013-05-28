#include <stdlib.h>
#include <iostream>
#include "grid.h"

int main(int argc, char* argv[])
{

std::cout << argv << std::endl;
  int L = atoi(argv[1]);
  
  std::cout << L << std::endl;
  
  Grid grid(L);
  
  return 0;
  
}