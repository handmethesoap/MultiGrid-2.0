#include "grid.h"
#include <iostream>
  
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
  _sigma = 0;
  _u = new double[_length_u];
  _f = new double[_length_f];

  for(int i = 0; i < _length_u; ++i)
  {
    _u[i] =  0;
  }

  for(int i = 0; i < _length_f; ++i)
  {
    _f[i] =  0;
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
		_u[i] =  0;
	}

	for(int i = 0; i < _length_f; ++i)
	{
		_f[i] =  0;
	}
}

Grid:: ~Grid()
{
  delete[] _u;
  delete[] _f;
}

void Grid:: set_u(int level, int x, int y, double value)
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

	_u[level_offset + (y)*dimension + x] = value;
	
}

double Grid:: get_u(int level, int x, int y)
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

	return _u[level_offset + (y)*dimension + x];
}

void Grid:: set_f(int level, int x, int y, double value)
{
	int dimension = (1 << level) - 1;
	int level_offset = 0;

	for(int i = 1; i < level; ++i)
  {
    level_offset += ((1 << i) - 1)*((1 << i) - 1);
  }

	if( x >= dimension )
	{
		std::cout << "max x dimension of " << dimension << " reached or execeed by attempt to access " << x << std::endl;
	}
	if( y >= dimension )
	{
		std::cout << "max y dimension of " << dimension << " reached or execeed by attempt to access " << y << std::endl;
	}

	_f[level_offset + (y)*dimension + x] = value;
	
}

double Grid:: get_f(int level, int x, int y)
{
	int dimension = (1 << level) - 1;
	int level_offset = 0;

	for(int i = 1; i < level; ++i)
  {
    level_offset += ((1 << i) - 1)*((1 << i) - 1);
  }

	if( x >= dimension )
	{
		std::cout << "max x dimension of " << dimension << " reached or execeed by attempt to access " << x << std::endl;
	}
	if( y >= dimension )
	{
		std::cout << "max y dimension of " << dimension << " reached or execeed by attempt to access " << y << std::endl;
	}

	return _f[level_offset + (y)*dimension + x];
}

void Grid:: initialise_u_boundary(double(* u_initialiser)(double, double))
{
  int dimension = (1 << _levels);
  for( int xy = 0; xy <= dimension; ++xy )
  {
    set_u(_levels, xy, 0, u_initialiser(xy/double(dimension), 0));
    set_u(_levels, xy, dimension, u_initialiser(xy/double(dimension), 1));
    set_u(_levels, 0, xy, u_initialiser(0, xy/double(dimension)));
    set_u(_levels, dimension, xy, u_initialiser(1, xy/double(dimension)));
  }
}

void Grid:: initialise_u(double(* u_initialiser)(double, double))
{
	int dimension = (1 << _levels);
	for( int y = 1; y < dimension; ++y )
	{
		for( int x = 1; x < dimension; ++x)
		{
			set_u(_levels, x, y, u_initialiser(x/double(dimension), y/double(dimension)));
		}
	}

}

void Grid:: initialise_f(double(* u_initialiser)(double, double))
{
	int dimension = (1 << _levels);
	for( int y = 0; y < dimension - 1; ++y )
	{
		for( int x = 0; x < dimension - 1; ++x)
		{
			set_f(_levels, x, y, u_initialiser((x+1)/double(dimension), (y+1)/double(dimension)));
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
	for(int y = 0; y <= (1 << level); ++y)
	{
		for( int x = 0; x <= (1 << level); ++x)
		{
			std::cout << get_u(level, x, y) << " ";
		}
		std::cout << std::endl;
	}
}

void Grid:: print_f(int level)
{
	for(int y = 0; y < (1 << level) - 1; ++y)
	{
		for( int x = 0; x < (1 << level) - 1; ++x)
		{
			std::cout << get_f(level, x, y) << " ";
		}
		std::cout << std::endl;
	}
}

void Grid:: rb_gauss_seidel_relaxation(int level)
{
  int dimension = (1 << level);
  double h2 = 1.0/double(dimension);
  int start;
  double value;
  double mult = 1.0/(4.0 + h2*_sigma);
  for( int y = 1; y < dimension; ++y )
  {
    start = y%2 + 1;
    for( int x = start; x < dimension; x += 2)
    {
      value = mult*(get_u(level, x+1, y) + get_u(level, x-1, y) + get_u(level, x, y+1) + get_u(level, x, y-1) + h2*get_f(level, x-1, y-1));
      set_u(level, x, y, value);
    }
  }

  for( int y = 1; y < dimension; ++y )
  {
    start = (y+1)%2 + 1;
    for( int x = start; x < dimension; x += 2)
    {
      value = mult*(get_u(level, x+1, y) + get_u(level, x-1, y) + get_u(level, x, y+1) + get_u(level, x, y-1) + h2*get_f(level, x-1, y-1));
      set_u(level, x, y, value);
    }
  }
}


void Grid:: fw_restrict(int level)
{
  double *residual = new double [((1<< level) + 1)*((1<< level) + 1)];
  
  
  delete residual;
}


void Grid:: calc_residual(double *residual)
{
  //set boundaries to zero
  
}










