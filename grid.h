#ifndef GRID_H
#define GRID_H

#include <iostream>

class Grid{

private:
  double* _u;
  double* _f;
  double* _r;
  unsigned int _levels, num_threads;
  unsigned int* _offsets_u;
  unsigned int* _offsets_f;

  unsigned int get_index(unsigned int level, unsigned int x, unsigned int y);
  unsigned int get_f_index(unsigned int level, unsigned int x, unsigned int y);

  unsigned int v1,v2;
public:
  Grid(unsigned int levels, unsigned int num_threads);
  ~Grid();

  friend std::ostream& operator<< (std::ostream &out, Grid &outputGrid);

  void rb_gauss_seidel_relaxation(unsigned int level,unsigned int times);
  void initialise_u_boundary(double(* u_initialiser)(double, double));
  void initialise_u(double(* u_initialiser)(double, double));
  void initialise_f(double(* u_initialiser)(double, double));

  void print(unsigned int level);
  void print_f(unsigned int level);

  void fw_restrict(unsigned int level);
  void calc_residual(unsigned int level,double &residual_norm);
  void calc_residual(unsigned int level);
  void interpolate(unsigned int level);
  void interpolate(unsigned int level, double(* boundary_initialiser)(double, double));

  void fmg_solve(unsigned int v_cycles[], unsigned int v1, unsigned int v2, double(* boundary_initialiser)(double, double));
  unsigned int Vcycle(unsigned int level, unsigned int v_cycles[], unsigned int v1, unsigned int v2);
  unsigned int Wcycle(unsigned int level, unsigned int w_cycles, unsigned int v1, unsigned int v2);

  double exact_L2_norm(unsigned int level, double(* exact_solution)(double, double));
  double L2_norm(unsigned int level);

};


#endif
