#include <stdlib.h>
#include <iostream>
#include "grid.h"
#include "initialisers.h"

int main(int argc, char* argv[])
{

  int L = atoi(argv[1]);
  double residual = 0.0;
  std::cout << L << std::endl;
  
  Grid A(L);
  A.initialise_u_boundary(constant);
	A.initialise_u_boundary(constant2);
  A.print(3);
	for(int i = 0; i < 100; ++i)
	{
		residual = A.rb_gauss_seidel_relaxation(3);
		std::cout << "residual = " << residual << std::endl;
	}
  A.print(3);
  A.print_f(3);
  
  
  return 0;
  
}