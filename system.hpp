#ifndef _system_hpp_
#define _system_hpp_
/*******************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.hpp"

/* функци€ индекса корневой €чейки */
int index_root(int i, int k, int l)
{
	return i*root_nz*root_ny+k*root_nz+l;
}

/* структура потока через границу */
struct flux
{
	double flux_rho, flux_rhovx, flux_rhovy, flux_rhovz, flux_rhoent;
};

/* структура €чейки вложенной сетки */
struct nested_cell	
{
	double rho,						// физические переменные
		   velx, vely, velz, 
		   entropy,
		   press, 
		   rhoent, rhovx, rhovy, rhovz;		// консервативные переменные
	
	// потоки через границу
	flux   x_plus, x_minus, y_plus, y_minus, z_plus, z_minus;
};

/* структура €чейки корневой сеткой */
struct root_cell		
{
	double rho,			// функци€ плотности - параметр мельчени€ сетки и потенциала
		   phi,			// гравитационный потенциал
		   phi_real,	// функции дл€ преобразовани€ ‘урье
		   phi_image;

	int nested_N;		// размер вложенной сетки
	double nested_h;	// шаг вложенной сетки [nested_h = root_h/nested_N]
	double *gravity;	// узлова€ вложенна€ сетка дл€ гравитации (nested_N+1)^3
	double *grav_new;	// узлова€ вложенна€ сетка дл€ итерацинного решател€
	nested_cell *mesh;	// вложенна€ сетка дл€ уравнений гидродинамики (nested_N)^3
	
	int index_nested_gravity_mesh(int i, int k, int l)	// индекс €чейки потенциала
	{
		return i*(nested_N+1)*(nested_N+1)+k*(nested_N+1)+l;
	}

	int index_nested_mesh(int i, int k, int l)	// индекс €чейки вложенной сетки
	{
		return i*nested_N*nested_N+k*nested_N+l;
	}

	void make_nested_mesh(int new_N)	// создание вложенной сетки
	{
		if(new_N > MAX_NESTED_N) new_N = MAX_NESTED_N;
		nested_N = new_N;
		nested_h = root_h/nested_N;
		mesh = new nested_cell[nested_N*nested_N*nested_N];
		gravity = new double[(nested_N+1)*(nested_N+1)*(nested_N+1)];
		grav_new = new double[(nested_N+1)*(nested_N+1)*(nested_N+1)];
	}

	void reconstruct_nested_mesh(int new_N)	// реконструкци€ вложенной сетки
	{
		return;		/*** функци€-заглушка ***/
	}

	void destroy_nested_mesh()	// удаление вложенной сетки
	{
		nested_N = 0;
		nested_h = 0.0;
		delete mesh;
		delete gravity;
		delete grav_new;
	}
};
/*******************************************************************/
#endif