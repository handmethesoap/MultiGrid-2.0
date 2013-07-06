#include "grid.h"
#include <iostream>
#include <fstream>
#include <cmath> 
  
Grid:: Grid(int levels)
{
  _offsets_u = new int[levels];
  _offsets_f = new int[levels]; 
  int temp = 0;
  int u_length_sum = 0;
  int f_length_sum = 0;

  for(int i = 1; i <= levels; ++i)
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
}

int Grid:: get_index(int level, int x, int y)
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

void Grid:: rb_gauss_seidel_relaxation(int level)
{
  int dimension = (1 << level);
  double h2 = 1.0/(double(dimension)*double(dimension));
  int start;
  double value;
  //double residual = 0.0;
  
  //loop through the even grid points
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
  for( int y = 1; y < dimension; ++y )
  {
    start = (y+1)%2 + 1;
    for( int x = start; x < dimension; x += 2)
    {
      value = 1.0/(4.0)*(_u[get_index(level, x+1, y)] + _u[get_index(level, x-1, y)] + _u[get_index(level, x, y+1)] + _u[get_index(level, x, y-1)] + h2*_f[get_f_index(level, x, y)]);
      //residual += std::abs(value - _u[get_index(level, x, y)]);
      _u[get_index(level, x, y)] = value;
    }
  }

}

void Grid:: fw_restrict(int level)
{
  int size = (1<< (level-1));
  int upper_size = (1<<level);
  double *residual = calc_residual(level);
  
  //caluclate restricted values
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      _f[get_f_index(level-1, x, y)] = 0.0625*( residual[(2*y-1)*(upper_size+1) + 2*x-1] + residual[(2*y-1)*(upper_size+1) + 2*x+1] + residual[(2*y+1)*(upper_size+1) + 2*x-1] + residual[(2*y+1)*(upper_size+1) + 2*x+1] + 2.0*( residual[(2*y-1)*(upper_size+1) + 2*x] + residual[(2*y+1)*(upper_size+1) + 2*x] + residual[2*y*(upper_size+1) + 2*x-1] + residual[2*y*(upper_size+1) + 2*x+1]) + 4.0*residual[2*y*(upper_size+1) + 2*x]);
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
  delete residual;
}

double * Grid:: calc_residual(int level)
{
  int size = (1 << level);
  double *residual = new double [(size + 1)*(size + 1)];
  double h2 = double(size)*double(size);
  
  //set boundary to zero
  for(int i = 0; i <= size; ++i)
  {
    residual[i] = 0.0;
    residual[(size +1)*size + i] = 0.0;
    residual[(size+1)*i] = 0.0;
    residual[(size+1)*i + size] = 0.0;
  }
  
  //caluclate interior points
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      residual[x + y*(size+1)] = h2*(_u[get_index(level, x + 1, y)] + _u[get_index(level, x - 1, y)] + _u[get_index(level, x, y + 1)] + _u[get_index(level, x, y - 1)] - 4.0*_u[get_index(level, x, y)]) + _f[get_f_index(level, x, y)];
    }
  }
  
  return residual;
  
}

void Grid:: interpolate(int level)
{
  int size = (1 << (level));
  int upper_size = (1 << (level + 1));
  double *interpol = new double [(upper_size + 1)*(upper_size + 1)];
  
  //set boundaries to zero
  for(int i = 0; i <= upper_size; ++i)
  {
    interpol[i] = 0.0;
    interpol[(upper_size)*(upper_size+1) + i] = 0.0;
    interpol[(upper_size+1)*i] = 0.0;
    interpol[(upper_size+1)*i + upper_size] = 0.0;
  }
  
  //copy directly values
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      interpol[2*y*(upper_size+1) + 2*x] = _u[get_index(level,x,y)];
    }
  }
  
  //horizontally interpolate values
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size + 1; ++x)
    {
      interpol[2*y*(upper_size+1) + 2*x - 1] = (interpol[2*y*(upper_size+1) + 2*x - 2] + interpol[2*y*(upper_size+1) + 2*x])*0.5;
    }
  }
  
  //vertically interpolate values
  for(int y = 1; y < size + 1; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      interpol[(2*y-1)*(upper_size+1) + 2*x] = (interpol[(2*y-2)*(upper_size+1) + 2*x] + interpol[2*y*(upper_size+1) + 2*x])*0.5;
    }
  }
 
  //diagonally interpolate values
  for(int y = 1; y < size + 1; ++y)
  {
    for(int x = 1; x < size + 1; ++x)
    {
      interpol[(2*y-1)*(upper_size+1) + 2*x-1] = (interpol[(2*y-2)*(upper_size+1) + 2*x-2] + interpol[(2*y-2)*(upper_size+1) + 2*x]+ interpol[(2*y)*(upper_size+1) + 2*x-2]+ interpol[(2*y)*(upper_size+1) + 2*x])*0.25;
    }
  }
 
  //add to gird
  for(int y = 0; y <= upper_size; ++y)
  {
    for(int x = 0; x <= upper_size; ++x)
    {
      _u[get_index(level+1, x,y)] += interpol[y*(upper_size+1) + x];
    }
  }
  
  delete interpol;
}

void Grid:: interpolate(int level, double(* boundary_initialiser)(double, double))
{
  int size = (1 << (level));
  int upper_size = (1 << (level + 1));
  
  //set boundaries
  for( int xy = 0; xy <= upper_size; ++xy )
  {
    _u[get_index(level+1, xy, 0)] =  boundary_initialiser(xy/double(upper_size), 0.0);
    _u[get_index(level+1, xy, upper_size)] = boundary_initialiser(xy/double(upper_size), 1.0);
    _u[get_index(level+1, 0, xy)] = boundary_initialiser(0.0, xy/double(upper_size));
    _u[get_index(level+1, upper_size, xy)] = boundary_initialiser(1.0, xy/double(upper_size));
  }
  
  //copy directly values
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      _u[get_index(level+1,2*x, 2*y)]= _u[get_index(level,x,y)];
    }
  }
  
  //horizontally interpolate values
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size + 1; ++x)
    {
      _u[get_index(level+1,2*x-1, 2*y)] = (_u[get_index(level+1,2*x-2, 2*y)] + _u[get_index(level+1,2*x, 2*y)])*0.5;
    }
  }
  
  //vertically interpolate values
  for(int y = 1; y < size + 1; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      _u[get_index(level+1,2*x, 2*y-1)] = (_u[get_index(level+1,2*x, 2*y-2)] + _u[get_index(level+1,2*x, 2*y)])*0.5;
    }
  }
 
  //diagonally interpolate values
  for(int y = 1; y < size + 1; ++y)
  {
    for(int x = 1; x < size + 1; ++x)
    {
      _u[get_index(level+1,2*x-1, 2*y-1)] = (_u[get_index(level+1,2*x-2, 2*y-2)] + _u[get_index(level+1,2*x, 2*y-2)]+ _u[get_index(level+1,2*x-2, 2*y)]+ _u[get_index(level+1,2*x, 2*y)])*0.25;
    }
  }
 
}

int Grid:: level_solve(int level, int v_cycles, int v1, int v2)
{
  int iterations = 0;
  double l2norm;
  double previousl2norm = 0.0;


  for(int k = 0; k < v_cycles; ++k)
  {
    //descend through the levels
    for( int i = level; i > 1; --i)
    {
      for( int j = 0; j < v1; ++j)
      {
	rb_gauss_seidel_relaxation(i);
      }
      fw_restrict(i);
    }

    //direct solution of lowest level
    rb_gauss_seidel_relaxation(1);

    //ascend through the levels
    for( int i = 1; i < level; ++i)
    {
      interpolate(i);
      for( int j = 0; j < v2; ++j)
      {
	rb_gauss_seidel_relaxation(i+1);
      }
    }
    ++iterations;
    #ifdef TEST
    l2norm = L2_norm(level);
    std::cout << "L2 norm of residual at level " << level << " = " << l2norm << std::endl;
    if(previousl2norm != 0.0)
    {
      std::cout << "convergence rate = " << (previousl2norm - l2norm)/previousl2norm << std::endl;
    }
    previousl2norm = l2norm;
    #endif
  }

  if(iterations == 1000)
  {
    return 0;
  }
  else
  {
    return iterations;
  }

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
  for( int y = 1; y < dimension; ++y )
  {
    for( int x = 1; x < dimension; ++x)
    {
      norm += (_u[get_index(level, x, y)] -  exact_solution(x/double(dimension), y/double(dimension)))*(_u[get_index(level, x, y)] -  exact_solution(x/double(dimension), y/double(dimension)));
    }
  }
  return sqrt(norm)/((dimension+1)*(dimension+1));
}

void Grid::fmg_solve(int v_cycles, double(* boundary_initialiser)(double, double))
{
  int cycles;
  double error;
  //set boundaries of lowest level
  int dimension = (1 << 1);
  for( int xy = 0; xy <= dimension; ++xy )
  {
    _u[get_index(1, xy, 0)] =  boundary_initialiser(xy/double(dimension), 0.0);
    _u[get_index(1, xy, dimension)] = boundary_initialiser(xy/double(dimension), 1.0);
    _u[get_index(1, 0, xy)] = boundary_initialiser(0.0, xy/double(dimension));
    _u[get_index(1, dimension, xy)] = boundary_initialiser(1.0, xy/double(dimension));
  }
  
  //solve lowest level
  rb_gauss_seidel_relaxation(1);

  for(int i = 1; i < _levels; ++i)
  {
    interpolate(i, boundary_initialiser);
    cycles = level_solve(i+1, v_cycles, 1, 1);
    
#ifdef TEST
    error = exact_L2_norm(i+1, boundary_initialiser);
    std::cout << "number of cycles at level " << i+1 << " = " << cycles << std::endl;
    std::cout << "L2 norm of error at level " << i+1 << " = " << error << std::endl << std::endl;
#endif
  }

}

double Grid::L2_norm(int level)
{
  double *residual = calc_residual(level);
  int size = (1 << level);
  double L2norm = 0.0;
  
  for(int y = 1; y < size; ++y)
  {
    for(int x = 1; x < size; ++x)
    {
      L2norm += residual[x + y*(size+1)] * residual[x + y*(size+1)];
    }
  }
  
  return sqrt(L2norm);
  
}












