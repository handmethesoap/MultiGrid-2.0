#ifndef GRID_H
#define GRID_H

#include <iostream>

class Grid{

private:
  double* _u;
  double* _f;
  double* _grid_pointers;
  int _levels;
  int* _offsets_u;
  int* _offsets_f;
  int* _grid_dimensions_u;
  int* _grid_dimensions_f;

  int get_index(int level, int x, int y);
  int get_f_index(int level, int x, int y);

  int v1,v2;
  double omega;
public:
  Grid(int levels);
  ~Grid();

  friend std::ostream& operator<< (std::ostream &out, Grid &outputGrid);

  void rb_gauss_seidel_relaxation(int level,int times);
  void initialise_u_boundary(double(* u_initialiser)(double, double));
  void initialise_u(double(* u_initialiser)(double, double));
  void initialise_f(double(* u_initialiser)(double, double));

  void print(int level);
  void print_f(int level);

  void fw_restrict(int level);
  double * calc_residual(int level,double &residual_norm);
  double * calc_residual(int level);
  void interpolate(int level);
  void interpolate(int level, double(* boundary_initialiser)(double, double));

  void fmg_solve(int v_cycles, int v1, int v2, double omega, double(* boundary_initialiser)(double, double));
  int Vcycle(int level, int v_cycles, int v1, int v2);
  int Wcycle(int level, int w_cycles, int v1, int v2);

  double exact_L2_norm(int level, double(* exact_solution)(double, double));
  double L2_norm(int level);

};


#endif
