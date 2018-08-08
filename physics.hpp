#ifndef _physics_hpp_
#define _physics_hpp_
/*******************************************************************/
#include <math.h>
#include "system.hpp"
#include "global.hpp"

#define ADIABATIC_INDEX_GAMMA (5.0/3.0)		// показатель адиабаты
#define MIN_RHO 0.0001						// минимальна€ плотность
#define KELVINS 0.1							// безразмерна€ температура

// функци€ начального профил€ плотности
double density_profile(double radius)
{
	if(radius < 1.0) 
		return 2.0*radius*radius*radius - 3.0*radius*radius + 1.0 + MIN_RHO;
	return MIN_RHO;
}

// функци€ начального профил€ давлени€
double pressure_profile(double radius)
{
	if(radius < 1.0) 
		return KELVINS*(2.0*radius*radius*radius - 3.0*radius*radius + 1.0 + MIN_RHO);
	return KELVINS*MIN_RHO;
}

// функци€ восстановлени€ энтропии
double get_entropy(double density, double pressure)
{
	return pressure/pow(density,ADIABATIC_INDEX_GAMMA);
}

// функци€ уравнени€ состо€ни€
double equation_of_state(double density, double entropy)
{
	return entropy * pow(density,ADIABATIC_INDEX_GAMMA);
}

// функци€ восстановлени€ физических переменных из консервативных
void Primitive_Variable()
{
	int i, k, l, i_nst, k_nst, l_nst, nx_nst, ny_nst, nz_nst;
	double dens, velx, vely, velz, press, rvelx, rvely, rvelz, entr, rhoent;
	double density, mass, h_nested;
	
	// восстановление физических переменных
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				nx_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				ny_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				nz_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				
				// работа с вложенной сеткой
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
						for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
						{
							// чтение плотности, момента импульса и плотности энтропии
							dens = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rho; 

							rvelx = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhovx;
								
							rvely = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhovy;
							
							rvelz = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhovz;

							rhoent = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhoent;
								
							velx = rvelx/dens;	// вычисление скорости
							vely = rvely/dens;
							velz = rvelz/dens;
							entr = rhoent/dens;	// вычисление энтропии

							// вычисление давлени€ по уравнению состо€ни€
							press = equation_of_state(dens,entr);		

							// запись скорости и давлени€
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].velx = velx;

							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].vely = vely;

							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].velz = velz;

							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].press = press;

							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].entropy = entr;
						}
			}

	// вычисление плотности в корневой €чейке
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				nx_nst   = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				ny_nst   = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				nz_nst   = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				h_nested = Root_Mesh_Hydro[index_root(i,k,l)].nested_h;
				
				// работа с вложенной сеткой
				mass = 0.0;
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
						for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
						{
							// чтение плотности
							density = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rho;
							
							// вычисление массы €чейки вложенной сетки
							mass += density*h_nested*h_nested*h_nested;
						}

				// вычисление и запись плотности корневой €чейки
				density = mass/root_h/root_h/root_h;
				Root_Mesh_Hydro[index_root(i,k,l)].rho = density;
			}
}
/*******************************************************************/
#endif