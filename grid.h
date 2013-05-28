#ifndef GRID_H
#define GRID_H

class Grid{

private:
  double* _u;
  double* _f;
	double _sigma;
	int _levels;
	int _length_u;
	int _length_f;

  void set_u(int level, int x, int y, double value);
	double get_u(int level, int x, int y);
	void set_f(int level, int x, int y, double value);
	double get_f(int level, int x, int y);

public:
  Grid(int levels);
  Grid(int levels, double sigma);
  ~Grid();
  
  void rb_gauss_seidel_relaxation(int level);
  void initialise_u_boundary(double(* u_initialiser)(double, double));
  void initialise_u(double(* u_initialiser)(double, double));
  void initialise_f(double(* u_initialiser)(double, double));
  void print(int level);
  void print_f(int level);
  void print_all(void);
  void print_all_f(void);
  void fw_restrict(int level);
  void calc_residual(double *residual);
  
  
};


#endif