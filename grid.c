#include "grid.h"
#include <iostream>
#include <cmath> 
  
Grid:: Grid(int levels)
{
  _length_u = 0;
  _length_f = 0; 
  int temp = 0;
  
  for(int i = 1; i <= levels; ++i)
  {
    temp = (1 << i);
    _length_u += (temp + 1)*(temp + 1);
    _length_f += (temp - 1)*(temp - 1);
  }
  
  _levels = levels;
  _sigma = 0.0;
  _u = new double[_length_u];
  _f = new double[_length_f];

  for(int i = 0; i < _length_u; ++i)
  {
    _u[i] =  0.0;
  }

  for(int i = 0; i < _length_f; ++i)
  {
    _f[i] =  0.0;
  }
}

Grid:: Grid(int levels, double sigma)
{
  _length_u = 0;
  _length_f = 0; 
  int temp = 0;
  
  for(int i = 1; i <= levels; ++i)
  {
    temp = (1 << i);
    _length_u += (temp + 1)*(temp + 1);
    _length_f += (temp - 1)*(temp - 1);
  }
  
	_levels = levels;
	_sigma = sigma;
  _u = new double[_length_u];
  _f = new double[_length_f];

	for(int i = 0; i < _length_u; ++i)
	{
		_u[i] =  0.0;
	}

	for(int i = 0; i < _length_f; ++i)
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

	for(int i = 1; i < level; ++i)
	{
	  level_offset += ((1 << i) + 1)*((1 << i) + 1);
	}

	if( x >= dimension )
	{
		std::cout << "max x dimension of " << dimension << " reached or execeed by attempt to access " << x << std::endl;
	}
	if( y >= dimension )
	{
		std::cout << "max y dimension of " << dimension << " reached or execeed by attempt to access " << y << std::endl;
	}

	return level_offset + (y)*dimension + x;
}

int Grid:: get_f_index(int level, int x, int y)
{
	int dimension = (1 << level) - 1;
	int level_offset = 0;

	for(int i = 1; i < level; ++i)
  {
    level_offset += ((1 << i) - 1)*((1 << i) - 1);
  }

	if( x-1 >= dimension )
	{
		std::cout << "max x dimension of " << dimension << " reached or execeed by attempt to access " << x << std::endl;
	}
	if( y-1 >= dimension )
	{
		std::cout << "max y dimension of " << dimension << " reached or execeed by attempt to access " << y << std::endl;
	}

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

void Grid:: print_all(void)
{
	for(int i = 0; i < _length_u; ++i)
	{
		std::cout << "_u[" << i+1 << "] = " << _u[i] << std::endl;
	}
	std::cout << "end" << std::endl;
}

void Grid:: print_all_f(void)
{
	for(int i = 0; i < _length_f; ++i)
	{
		std::cout << "_f[" << i+1 << "] = " << _f[i] << std::endl;
	}
	std::cout << "end" << std::endl;
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

double Grid:: rb_gauss_seidel_relaxation(int level)
{
  int dimension = (1 << level);
  double h2 = 1.0/(double(dimension)*double(dimension));
  int start;
  double value;
  double mult = 1.0/(4.0 + h2*_sigma);
  double residual = 0.0;
  
  for( int y = 1; y < dimension; ++y )
  {
    start = y%2 + 1;
    for( int x = start; x < dimension; x += 2)
    {
      value = mult*(_u[get_index(level, x+1, y)] + _u[get_index(level, x-1, y)] + _u[get_index(level, x, y+1)] + _u[get_index(level, x, y-1)] + h2*_f[get_f_index(level, x, y)]);
      residual += std::abs(value - _u[get_index(level, x, y)]);
      _u[get_index(level, x, y)] =  value;
			
    }
  }

  for( int y = 1; y < dimension; ++y )
  {
    start = (y+1)%2 + 1;
    for( int x = start; x < dimension; x += 2)
    {
      value = mult*(_u[get_index(level, x+1, y)] + _u[get_index(level, x-1, y)] + _u[get_index(level, x, y+1)] + _u[get_index(level, x, y-1)] + h2*_f[get_f_index(level, x, y)]);
      residual += std::abs(value - _u[get_index(level, x, y)]);
      //std::cout << std::abs(value - _u[get_index(level, x, y)]) << std::endl;
      _u[get_index(level, x, y)] = value;
    }
  }
  
	return residual*h2;
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
//   int x=1;
//   int y=1;
//   std::cout << "0.625*("<<residual[(2*y-1)*(upper_size+1) + 2*x-1]<< "+" << residual[(2*y-1)*(upper_size+1) + 2*x+1] << "+" << residual[(2*y+1)*(upper_size+1) + 2*x-1] << "+" << residual[(2*y+1)*(upper_size+1) + 2*x+1] << "+" << "2.0*(" << residual[(2*y-1)*(upper_size+1) + 2*x] << "+" << residual[(2*y+1)*(upper_size+1) + 2*x] <<"+"<<  residual[2*y*(upper_size+1) + 2*x-1] <<"+"<< residual[2*y*(upper_size+1) + 2*x+1]<<") +" << "4.0*" << residual[2*y*(upper_size+1) + 2*x] << ") = "<< 0.0625*( residual[(2*y-1)*(upper_size+1) + 2*x-1] + residual[(2*y-1)*(upper_size+1) + 2*x+1] + residual[(2*y+1)*(upper_size+1) + 2*x-1] + residual[(2*y+1)*(upper_size+1) + 2*x+1] + 2.0*( residual[(2*y-1)*(upper_size+1) + 2*x] + residual[(2*y+1)*(upper_size+1) + 2*x] + residual[2*y*(upper_size+1) + 2*x-1] + residual[2*y*(upper_size+1) + 2*x+1]) + 4.0*residual[2*y*(upper_size+1) + 2*x]) << std::endl;
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
//   int x=2;
//   int y=1;
//   std::cout << h2 << "*(" << _u[get_index(level, x + 1, y)] << "+" << _u[get_index(level, x - 1, y)] << "+" << _u[get_index(level, x, y + 1)] << "+" << _u[get_index(level, x, y - 1)] <<"- 4*" << 4.0*_u[get_index(level, x, y)] << ") +" << _f[get_f_index(level, x, y)] << "=" << h2*(_u[get_index(level, x + 1, y)] + _u[get_index(level, x - 1, y)] + _u[get_index(level, x, y + 1)] + _u[get_index(level, x, y - 1)] - 4.0*_u[get_index(level, x, y)]) + 0.0<< std::endl;
  
  //test print residual
//   for(int y = 0; y <= size; ++y)
//   {
//     for(int x = 0; x <= size; ++x)
//     {
//       std::cout << residual[x + y*(size+1)] << " ";
//     }
//     std::cout << std::endl;
//   }
  
  return residual;
  
}

void Grid:: interpolate(int level)
{
  int size = (1 << (level));
  int upper_size = (1 << (level + 1));
  double *interpol = new double [(upper_size + 1)*(upper_size + 1)];
  
  //temporarily set all values to 3
  for(int y = 0; y <= upper_size; ++y)
  {
    for(int x = 0; x <= upper_size; ++x)
    {
      interpol[y*(upper_size+1) + x] = 3;
    }
  }
  
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
  
//   //print interpolated grid
//   for(int y = 0; y <= upper_size; ++y)
//   {
//     for(int x = 0; x <= upper_size; ++x)
//     {
//       std::cout << interpol[y*(upper_size + 1) + x] << " ";
//     }
//     std::cout << std::endl;
//   }
 
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








