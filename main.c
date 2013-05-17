#include <stdlib.h>
#include <iostream>
#include "grid.h"

int main(int argc, char* argv[])
{
  int L = atoi(argv[1]);
  
  std::cout << L << std::endl;
  
  Grid grid(L);
  
  return 0;
  
}