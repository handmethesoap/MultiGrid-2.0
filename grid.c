#include "grid.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;
Grid:: Grid(int levels, int _num_threads)
{
  _offsets_u = new int[levels];
  _offsets_f = new int[levels];
  int temp = 0;
  int u_length_sum = 0;
  int f_length_sum = 0;

  for (int i = 1; i <= levels; ++i)
  {
    _offsets_u[i-1] = u_length_sum;
    _offsets_f[i-1] = f_length_sum;
    temp = (1 << i);
    u_length_sum += (temp + 1)*(temp + 1);
    f_length_sum += (temp - 1)*(temp - 1);
  }

  _levels = levels;
  _u = new double[u_length_sum];
  _f = new double[f_length_sum];
  _r = new double[u_length_sum];
  num_threads = _num_threads;
  for(int i = 0; i < u_length_sum; ++i)
  {
    _u[i] =  0.0;
  }

  for(int i = 0; i < f_length_sum; ++i)
  {
    _f[i] =  0.0;
  }
}

Grid:: ~Grid()
{
  delete[] _u;
  delete[] _f;
  delete[] _r;
  delete[] _offsets_u;
  delete[] _offsets_f;
}

int Grid:: get_index(int level, int x,int y)
{
  int dimension = (1 << level) + 1;
  int level_offset = 0;

  level_offset = _offsets_u[level - 1];

#ifdef DEBUG
  if( x >= dimension )
  {
    std::cout << "max x dimension of " << dimension << " reached or execeed by attempt to access " << x << std::endl;
  }
  if( y >= dimension )
  {
    std::cout << "max y dimension of " << dimension << " reached or execeed by attempt to access " << y << std::endl;
  }
#endif

  return level_offset + (y)*dimension + x;
}

int Grid:: get_f_index(int level, int x, int y)
{
  int dimension = (1 << level) - 1;
  int level_offset = 0;

  level_offset = _offsets_f[level - 1];

#ifdef DEBUG
  if( x-1 >= dimension )
  {
	  std::cout << "max x dimension of " << dimension << " reached or execeed by attempt to access " << x << std::endl;
  }
  if( y-1 >= dimension )
  {
	  std::cout << "max y dimension of " << dimension << " reached or execeed by attempt to access " << y << std::endl;
  }
#endif

  return level_offset + (y - 1)*dimension + (x - 1);
}

void Grid:: initialise_u_boundary(double(* u_initialiser)(double, double))
{
  int dimension = (1 << _levels);
  for( int xy = 0; xy <= dimension; ++xy )
  {
    _u[get_index(_levels, xy, 0)] =  u_initialiser(xy/double(dimension), 0.0);
    _u[get_index(_levels, xy, dimension)] = u_initialiser(xy/double(dimension), 1.0);
    _u[get_index(_levels, 0, xy)] = u_initialiser(0.0, xy/double(dimension));
    _u[get_index(_levels, dimension, xy)] = u_initialiser(1.0, xy/double(dimension));
  }
}

void Grid:: initialise_u(double(* u_initialiser)(double, double))
{
  int dimension = (1 << _levels);
  for( int y = 1; y < dimension; ++y )
  {
    for( int x = 1; x < dimension; ++x)
    {
	    _u[get_index(_levels, x, y)] =  u_initialiser(x/double(dimension), y/double(dimension));
    }
  }

}

void Grid:: initialise_f(double(* u_initialiser)(double, double))
{
  int dimension = (1 << _levels);
  for( int y = 1; y < dimension; ++y )
  {
    for( int x = 1; x < dimension; ++x)
    {
	    _f[get_f_index(_levels, x, y)] = u_initialiser((x)/double(dimension), (y)/double(dimension));
    }
  }

}

void Grid:: print(int level)
{
  std::cout << std::endl;
  for(int y = 0; y <= (1 << level); ++y)
  {
    for( int x = 0; x <= (1 << level); ++x)
    {
      std::cout << _u[get_index(level, x, y)] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void Grid:: print_f(int level)
{
  std::cout << std::endl;
  for(int y = 1; y <= (1 << level) - 1; ++y)
  {
    for( int x = 1; x <= (1 << level)- 1; ++x)
    {
      std::cout << _f[get_f_index(level, x, y)] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void Grid:: rb_gauss_seidel_relaxation(int level,int times)
{
    int dimension = (1 << level);
    double h2 = 1.0/(double(dimension)*double(dimension));
    int start;
    double value;
    //double residual = 0.0;
    for (int iter = 0; iter<times; ++iter){
        //loop through the even grid points
    #ifdef USEOMP
        #pragma omp parallel for private(start,value) if(num_threads>1)
    #endif
        for( int y = 1; y < dimension; ++y )
        {
            start = y%2 + 1;
            for( int x = start; x < dimension; x += 2)
            {
                value = 1.0/(4.0)*(_u[get_index(level, x+1, y)] + _u[get_index(level, x-1, y)] + _u[get_index(level, x, y+1)] + _u[get_index(level, x, y-1)] + h2*_f[get_f_index(level, x, y)]);
                //residual += std::abs(value - _u[get_index(level, x, y)]);
                _u[get_index(level, x, y)] =  value;
            }
        }
        //loop through the odd grid points
    #ifdef USEOMP
        #pragma omp parallel for private(start,value) if(num_threads>1)
    #endif
        for( int y = 1; y < dimension; ++y )
        {
            start = (y+1)%2 + 1;
            for( int x = start; x < dimension; x += 2)
            {
                value = 1.0/(4.0)*(_u[get_index(level, x+1, y)] + _u[get_index(level, x-1, y)] + _u[get_index(level, x, y+1)] + _u[get_index(level, x, y-1)] + h2*_f[get_f_index(level, x, y)]);
                _u[get_index(level, x, y)] = value;
            }
        }
    }
}

void Grid:: fw_restrict(int level)
{
  int size = (1<< (level-1));
  calc_residual(level);
  //caluclate restricted values
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      _f[get_f_index(level-1, x, y)] = 0.0625*( _r[get_index(level, 2*x-1, 2*y-1)] +
                                               _r[get_index(level, 2*x + 1, 2*y -1)] +
                                               _r[get_index(level, 2*x-1, 2*y+1)] +
                                               _r[get_index(level, 2*x+1, 2*y+1)] +
                                               2.0*( _r[get_index(level, 2*x, 2*y-1)] +
                                                    _r[get_index(level, 2*x, 2*y+1)] +
                                                    _r[get_index(level, 2*x-1, 2*y)] +
                                                    _r[get_index(level, 2*x+1, 2*y)]) +
                                               4.0*_r[get_index(level, 2*x, 2*y)]);
    }
  }

  //write zeroes to the _u matrix for the level we just restircted to
  for(int y = 0; y <= size; ++y)
  {
    for(int x = 0; x <= size; ++x)
    {
      _u[get_index(level-1, x, y)] = 0.0;
    }
  }
}

void Grid:: calc_residual(int level)
{
  int size = (1 << level);
  double h2 = double(size)*double(size);

  //set boundary to zero
  for( int xy = 0; xy <= size; ++xy )
  {
    _r[get_index(level, xy, 0)] =  0.0;
    _r[get_index(level, xy, size)] = 0.0;
    _r[get_index(level, 0, xy)] = 0.0;
    _r[get_index(level, size, xy)] = 0.0;
  }

  //caluclate interior points
#ifdef USEOMP
  #pragma omp parallel for if(num_threads>1)
#endif
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      _r[get_index(level, x, y)] = h2*(_u[get_index(level, x + 1, y)] + _u[get_index(level, x - 1, y)] + _u[get_index(level, x, y + 1)] + _u[get_index(level, x, y - 1)] - 4.0*_u[get_index(level, x, y)]) + _f[get_f_index(level, x, y)];
    }
  }

}

void Grid::interpolate(int level){
  int size = (1 << (level));
  //upper row
#ifdef USEOMP
  #pragma omp parallel for if(num_threads>1)
#endif
  for(int y = 1; y < size; ++y){
    for(int x = 1; x < size; ++x){
      int xx = 2*x, yy = 2*y;
      _u[get_index(level+1, xx-1,yy+1)] += 0.25 * _u[get_index(level,x,y)];
      _u[get_index(level+1, xx  ,yy+1)] += 0.5  * _u[get_index(level,x,y)];
      _u[get_index(level+1, xx+1,yy+1)] += 0.25 * _u[get_index(level,x,y)];
      _u[get_index(level+1, xx-1,yy)]   += 0.5  * _u[get_index(level,x,y)];
      _u[get_index(level+1, xx,  yy)]   +=        _u[get_index(level,x,y)];
      _u[get_index(level+1, xx+1,yy)]   += 0.5  * _u[get_index(level,x,y)];
      _u[get_index(level+1, xx-1,yy-1)] += 0.25 * _u[get_index(level,x,y)];
      _u[get_index(level+1, xx  ,yy-1)] += 0.5  * _u[get_index(level,x,y)];
      _u[get_index(level+1, xx+1,yy-1)] += 0.25 * _u[get_index(level,x,y)];
    }
  }

}

void Grid:: interpolate(int level, double(* boundary_initialiser)(double, double))
{
  int size = (1 << (level));
  int upper_size = (1 << (level + 1));

//set boundaries
for( int xy = 0; xy <= upper_size; ++xy ){
    _u[get_index(level+1, xy, 0)] =  boundary_initialiser(xy/double(upper_size), 0.0);
    _u[get_index(level+1, xy, upper_size)] = boundary_initialiser(xy/double(upper_size), 1.0);
    _u[get_index(level+1, 0, xy)] = boundary_initialiser(0.0, xy/double(upper_size));
    _u[get_index(level+1, upper_size, xy)] = boundary_initialiser(1.0, xy/double(upper_size));
}


        //copy directly values
        for(int y = 1; y < size; ++y){
            for(int x = 1; x < size; ++x){
                _u[get_index(level+1,2*x, 2*y)]= _u[get_index(level,x,y)];
            }
        }

        //horizontally interpolate values
        for(int y = 1; y < size; ++y){
            for(int x = 1; x < size + 1; ++x){
                //_u[get_index(level+1,2*x-1, 2*y)] = (_u[get_index(level+1,2*x-2, 2*y)] + _u[get_index(level+1,2*x, 2*y)])*0.5;
                _u[get_index(level+1,2*x-1, 2*y)] = (_u[get_index(level,x-1, y)] + _u[get_index(level,x, y)])*0.5;
            }
        }

        //vertically interpolate values
        for(int y = 1; y < size + 1; ++y){
            for(int x = 1; x < size; ++x){
                //_u[get_index(level+1,2*x, 2*y-1)] = (_u[get_index(level+1,2*x, 2*y-2)] + _u[get_index(level+1,2*x, 2*y)])*0.5;
                _u[get_index(level+1,2*x, 2*y-1)] = (_u[get_index(level,x, y-1)] + _u[get_index(level,x, y)])*0.5;
            }
        }

        //diagonally interpolate values
        for(int y = 1; y < size + 1; ++y){
            for(int x = 1; x < size + 1; ++x){
                // _u[get_index(level+1,2*x-1, 2*y-1)] = (_u[get_index(level+1,2*x-2, 2*y-2)] + _u[get_index(level+1,2*x, 2*y-2)]+ _u[get_index(level+1,2*x-2, 2*y)]+ _u[get_index(level+1,2*x, 2*y)])*0.25;
                _u[get_index(level+1,2*x-1, 2*y-1)] = (_u[get_index(level,x-1, y-1)] + _u[get_index(level,x, y-1)]+ _u[get_index(level,x-1, y)]+ _u[get_index(level,x, y)])*0.25;
            }
        }

}


int Grid::Vcycle(int level, int v_cycles, int v1, int v2)
{
     int iteration = 0;
     //#ifdef TEST
	 double previousl2norm = 0.0;
     //#endif
    for(int k = 0; k < v_cycles; ++k)
    {
        //descend through the levels
        for( int i = level; i > 1; --i)
        {
          rb_gauss_seidel_relaxation(i,v1);
          fw_restrict(i);
        }

        //direct solution of lowest level
        rb_gauss_seidel_relaxation(1,1);

        //ascend through the levels
        for(int i = 1; i < level; ++i)
        {
            interpolate(i);
            rb_gauss_seidel_relaxation(i+1,v2);
        }

        #ifdef TEST

        double l2norm = L2_norm(level);
        std::cout << "L2 norm of residual at level " << level << " = " << l2norm << std::endl;

        if(previousl2norm != 0.0)
        {
            std::cout << "convergence rate = " << (previousl2norm - l2norm)/previousl2norm << std::endl;
        }
        previousl2norm = l2norm;
        #endif

        ++iteration;
    }
    return iteration;
}

std::ostream& operator<< (std::ostream &out, Grid &outputGrid)
{
  double n = ( 1<< outputGrid._levels);

  //output u matrix in format compatiable with gnuplot
  for( int y = 0; y < n; ++y)
  {
    for( int x = 0; x < n; ++x )
    {
      out << double(x)/(n-1.0) << " " <<  double(y)/(n-1.0) << " " << outputGrid._u[outputGrid.get_index(outputGrid._levels, x, y)] << std::endl;
    }
    out << std::endl;
  }

  return out;
}

double Grid:: exact_L2_norm(int level, double(* exact_solution)(double, double))
{
  int dimension = (1 << level);
  double norm = 0.0;
#ifdef USEOMP
  #pragma omp parallel for reduction(+:norm) if(num_threads>1)
#endif
  for( int y = 1; y < dimension; ++y )
  {
    for(int x = 1; x < dimension; ++x)
    {
      norm += (_u[get_index(level, x, y)] -  exact_solution(x/double(dimension), y/double(dimension)))*(_u[get_index(level, x, y)] -  exact_solution(x/double(dimension), y/double(dimension)));
    }
  }
  return sqrt(norm)/((dimension+1)*(dimension+1));
}

void Grid::fmg_solve(int v_cycles, int _v1,int _v2, double(* boundary_initialiser)(double, double))
{

  #ifdef USEOMP
    omp_set_num_threads(num_threads);
  #endif // USEOMP
  //set boundaries of lowest level
  int dimension = (1 << 1);
  v1 = _v1;
  v2 = _v2;
  for( int xy = 0; xy <= dimension; ++xy )
  {
    _u[get_index(1, xy, 0)] =  boundary_initialiser(xy/double(dimension), 0.0);
    _u[get_index(1, xy, dimension)] = boundary_initialiser(xy/double(dimension), 1.0);
    _u[get_index(1, 0, xy)] = boundary_initialiser(0.0, xy/double(dimension));
    _u[get_index(1, dimension, xy)] = boundary_initialiser(1.0, xy/double(dimension));
  }

  //solve lowest level
  rb_gauss_seidel_relaxation(1,1);

  for(int i = 1; i < _levels; ++i)
  {
    interpolate(i, boundary_initialiser);
    Vcycle(i+1, v_cycles, v1, v2);

#ifdef TEST
    double error = exact_L2_norm(i+1, boundary_initialiser);
    std::cout << "L2 norm of error at level " << i+1 << " = " << error << std::endl << std::endl;
#endif
  }

}


double Grid::L2_norm(int level)
{
  calc_residual(level);
  int size = (1 << level);
  double L2norm = 0.0;
#ifdef USEOMP
  #pragma omp parallel for reduction(+:L2norm) if(num_threads>1)
#endif
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      L2norm += _r[get_index(level, x, y)] * _r[get_index(level, x, y)];
    }
  }

  return sqrt(L2norm);

}












