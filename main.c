// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <sys/time.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include "grid.h"
#include "initialisers.h"

int main(int argc, char* argv[])
{
  
  //initialise timing variables
  timeval start;
  timeval end;
  gettimeofday(&start, NULL);
  
  //Initialse variables
  int operation_time = 0;
  
  int L = atoi(argv[1]);
  int it = 0;

  ////initialise grid
  Grid A(L);
  A.initialise_u_boundary(sinh_function);
  A.initialise_sigma(zero);


  //solve multigrid problem and output number of iterations
  it = A.multigrid_solve(1e-10, 1, 1);
  std::cout << "number of iterations = "<< it << std::endl;

  //end time measurement
  gettimeofday(&end, NULL);
  
  //output result to file in format compatiable with gnuplot
  std::ofstream outputfile;
  outputfile.open("solution.dat");
  outputfile << A;
  outputfile.close();

  //calculate and output operating time
  operation_time = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec; 
  std::cout << "Operating time (us): " << operation_time << std::endl;
  std::cout << "start time " << start.tv_sec << ":" << start.tv_usec << std::endl;
  std::cout << "end time " << end.tv_sec << ":" << end.tv_usec << std::endl;
  
  #ifndef DEBUG
  //log operating time
  time_t now = time(0);
  tm *ltm = localtime(&now);
  std::ofstream logfile;
  logfile.open("performance_log.txt",  std::ofstream::out | std::ofstream::app);
  logfile << ltm->tm_mday << "/" << 1 + ltm->tm_mon << "/" << 1900 + ltm->tm_year << " " << 1 + ltm->tm_hour << ":" << 1 + ltm->tm_min << ":" << 1 + ltm->tm_sec << " - " << operation_time << " us" << std::endl;
  logfile.close();
  #endif

  return 0;
}
