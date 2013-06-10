// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "grid.h"
#include "initialisers.h"

int main(int argc, char* argv[])
{
  int L = atoi(argv[1]);
  int it = 0;

  ////initialise grid
  Grid A(L);
  A.initialise_u_boundary(sin_function);
  A.initialise_f(sin_function_2);
  A.initialise_sigma(omega);


  //solve multigrid problem and output number of iterations
  it = A.multigrid_solve(1e-10, 1, 1);
  std::cout << "number of iterations = "<< it << std::endl;

  //output result to file in format compatiable with gnuplot
  std::ofstream outputfile;
  outputfile.open("solution.dat");
  outputfile << A;
  outputfile.close();


  return 0;
}
