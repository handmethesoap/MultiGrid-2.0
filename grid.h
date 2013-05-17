#ifndef GRID_H
#define GRID_H

class Grid{

private:
  double* _u;
  double* _f;
  
public:
  Grid(int levels);
  ~Grid();
  
  void set_u_boundary(double(* u_initialiser)(double, double));
  
  
};


#endif