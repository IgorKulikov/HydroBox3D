#ifndef _cfl_hpp_
#define _cfl_hpp_
/*******************************************************************/
#include "physics.hpp"
#include "global.hpp"

// Определение шага по времени по числу Куранта
void Tau_Computing(double timer, double tend)
{
	int i, k, l, i_nst, k_nst, l_nst, nx_nst, ny_nst, nz_nst;
	double dens, velx, vely, velz, vel, press; 
	double min_h = root_h, sound, actual_vel, actual_h, max_vel = 0.0;

	// цикл по корневой сетке
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				nx_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				ny_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				nz_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				actual_h = Root_Mesh_Hydro[index_root(i,k,l)].nested_h;
		
				// определение минимального шага по пространству
				if(actual_h < min_h) min_h = actual_h;
				
				// работа с вложенной сеткой
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
						for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
						{
							// чтение плотности, скорости и давления
							dens = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rho; 

							velx = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].velx; 

							vely = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].vely;

							velz = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].velz;

							press = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].press;
								
							// вычисление актуальной скорости
							vel = sqrt(velx*velx+vely*vely+velz*velz);
							sound = sqrt(ADIABATIC_INDEX_GAMMA*press/dens);
							actual_vel = vel + sound;

							// выделение максимальной скорости в расчетной области
							if(actual_vel > max_vel) max_vel = actual_vel;
						}
			}

	// Вычисление временного шага из условия Куранта
	tau = min_h * CFL / max_vel;
	if(timer + tau > tend) tau = tend - timer;
}
/*******************************************************************/
#endif