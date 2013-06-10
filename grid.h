#ifndef GRID_H
#define GRID_H

#include <iostream>

class Grid{

private:
  double* _u;
  double* _f;
	double* _sigma;
	int _levels;
	int _length_u;
	int _length_f;

	int get_index(int level, int x, int y);
	int get_f_index(int level, int x, int y);

public:
  Grid(int levels);
  ~Grid();

  friend std::ostream& operator<< (std::ostream &out, Grid &outputGrid);

  double rb_gauss_seidel_relaxation(int level);

  void initialise_u_boundary(double(* u_initialiser)(double, double));
  void initialise_u(double(* u_initialiser)(double, double));
  void initialise_f(double(* u_initialiser)(double, double));
	void initialise_sigma(double(* u_initialiser)(double, double));

  void print(int level);
  void print_f(int level);
  void print_all(void);
	void print_all_f(void);
	void print_sigma(int level);

  void fw_restrict(int level);
  double * calc_residual(int level);
  void interpolate(int level);

	int rb_solve(double precision);
	int two_level_solve(double precision, int v1, int v2);
	int multigrid_solve(double precision, int v1, int v2);

	double L2_norm(double(* exact_solution)(double, double));

	

  
  
};


#endif