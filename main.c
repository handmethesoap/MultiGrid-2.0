#include <stdlib.h>
#include <iostream>
#include "grid.h"
#include "initialisers.h"

int main(int argc, char* argv[])
{

  int L = atoi(argv[1]);
  
  std::cout << L << std::endl;
  
  Grid A(L);
  A.initialise_u_boundary(constant);
  A.print(3);
  A.rb_gauss_seidel_relaxation(3);
  A.print(3);
  A.print_f(3);
  
  
  return 0;
  
}