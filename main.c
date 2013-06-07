#include <stdlib.h>
#include <iostream>
#include "grid.h"
#include "initialisers.h"

int main(int argc, char* argv[])
{

  int L = atoi(argv[1]);
  double residual = 0.0;
  int it = 0;
  std::cout << L << std::endl;
  
  Grid A(L);
  A.initialise_u_boundary(constant);
  //A.initialise_u_boundary(constant2);
//   for(int i = 0; i <= 5; ++i)
//   {
//     std::cout << "----------------------ITERATION NUMBER " << i << std::endl;
//     A.print(3);
//     A.rb_gauss_seidel_relaxation(3);
//     A.print(3);
//     A.fw_restrict(3);
//     A.print_f(2);
//     A.rb_gauss_seidel_relaxation(2);
//     A.print(2);
//     A.interpolate(2);
//     A.print(3);
//     
//   }
  do
  {
    residual = A.rb_gauss_seidel_relaxation(3);
    std::cout << "residual = " << residual << std::endl;
    A.print(3);
    ++it;
    A.fw_restrict(3);
    residual = A.rb_gauss_seidel_relaxation(2);
    A.fw_restrict(2);
    residual = A.rb_gauss_seidel_relaxation(1);
    A.interpolate(1);
    ++it;
    residual = A.rb_gauss_seidel_relaxation(2);
    A.interpolate(2);
    A.print(3);
    ++it;
    
  }while(residual > 1e-7);
  std::cout << "number of iterations = "<< it << std::endl;
  A.print(3);

  
  return 0;
  
}