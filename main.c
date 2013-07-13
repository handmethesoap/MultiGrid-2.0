// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <sys/time.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include "grid.h"
#include "initialisers.h"

using namespace std;

int main(int argc, char* argv[])
{
  int L = 10;
  unsigned num_threads = 1;
  stringstream ss("");
  if (argc>3){
    cerr<<"Too many arguments. Syntax: ./multigrid num_levels num_threads"<<endl;
    return -1;
  }
  if (argc>2){
    ss<<argv[1]<<" "<<argv[2];
    ss>>L>>num_threads;
    if (ss.fail()){
        cerr<<"Error reading arguments.Syntax: ./multigrid num_levels num_threads"<<endl;
        return -1;
    }
  }
  else if (argc>1){
    ss<<argv[1];
    ss>>L;
    if (ss.fail()){
        cerr<<"Error reading arguments.Syntax: ./multigrid num_levels num_threads"<<endl;
        return -1;
    }
  }
  cout<<"Running multigrid with"<<endl;
  cout<<"Number of levels L = "<<L<<endl;
  cout<<"Number of threads = "<<num_threads<<endl;
  //initialise timing variables
  timeval start;
  timeval end;
  gettimeofday(&start, NULL);

  //Initialse variables
  int operation_time = 0;


  ////initialise grid
  Grid A(L,num_threads);
  A.initialise_u_boundary(sinh_function);

  //int vcount[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  A.fmg_solve(2,5,5,sinh_function);

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

  return 0;
}
