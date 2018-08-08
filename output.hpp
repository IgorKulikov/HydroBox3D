#ifndef _output_hpp_
#define _output_hpp_
/*******************************************************************/
#include <stdio.h>
#include "global.hpp"
#include "physics.hpp"

// ‘ункци€ сохранени€ величин энергий
void LogWrite(double timer)
{
	int i, k, l, i_nst, k_nst, l_nst, nx_nst, ny_nst, nz_nst;
	double h_nst, volume, density, velocity, velx, vely, velz, potential, pressure;
	double mass, energy_kinetic, energy_gravity, energy_internal;

	// инициализаци€ интегралов
	mass = 0.0;
	energy_kinetic  = 0.0; 
	energy_gravity  = 0.0;
	energy_internal = 0.0;

	// цикл по корневой сетке
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				nx_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				ny_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				nz_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				h_nst  = Root_Mesh_Hydro[index_root(i,k,l)].nested_h;
				volume = h_nst * h_nst * h_nst;

				// цикл по вложенной сетке
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
						for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
						{
							// чтение параметров
							density = Root_Mesh_Hydro[index_root(i,k,l)].
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
							velocity = sqrt(velx*velx + vely*vely + velz*velz);
							
							potential = 0.125 * ( 
							 Root_Mesh_Hydro[index_root(i,k,l)].
							  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
							  index_nested_gravity_mesh(i_nst+1,k_nst+1,l_nst+1)] + 
							 Root_Mesh_Hydro[index_root(i,k,l)].
							  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
							  index_nested_gravity_mesh(i_nst+1,k_nst+1,l_nst  )] +
							 Root_Mesh_Hydro[index_root(i,k,l)].
							  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
							  index_nested_gravity_mesh(i_nst+1,k_nst  ,l_nst+1)] +
							 Root_Mesh_Hydro[index_root(i,k,l)].
							  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
							  index_nested_gravity_mesh(i_nst+1,k_nst  ,l_nst  )] +
							 Root_Mesh_Hydro[index_root(i,k,l)].
							  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
							  index_nested_gravity_mesh(i_nst  ,k_nst+1,l_nst+1)] + 
							 Root_Mesh_Hydro[index_root(i,k,l)].
							  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
							  index_nested_gravity_mesh(i_nst  ,k_nst+1,l_nst  )] +
							 Root_Mesh_Hydro[index_root(i,k,l)].
							  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
							  index_nested_gravity_mesh(i_nst  ,k_nst  ,l_nst+1)] +
							 Root_Mesh_Hydro[index_root(i,k,l)].
							  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
							  index_nested_gravity_mesh(i_nst  ,k_nst  ,l_nst  )]);

							pressure = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].press;

							// добавление в интеграл значений энергий и массы
							mass += density*volume;
							energy_kinetic  += density*velocity*velocity/2.0*volume; 
							energy_gravity  += density*potential/2.0*volume;
							energy_internal += 
								pressure/(ADIABATIC_INDEX_GAMMA-1.0)*volume;
						}
			}

	// сохранение результатов в файл
	fprintf(fout,"%lf %lf %lf %lf %lf %lf\n",
		timer, mass, energy_kinetic, energy_gravity, energy_internal,
		energy_kinetic + energy_gravity + energy_internal);
	fflush(fout);
}

// ‘ункци€ сохранени€ результата
void Save(double timer)
{
	int MAX_N = MAX_NESTED_N * N, l_nested_nz;
	int i, k, i_root, k_root, l_root, i_nested, k_nested, l_nested; 
	double *column_density; 
	double mass_column, density, total_density, h_nested, h_print, x, y;	
	FILE *fout;
	char FileName[128];
	
	column_density = new double[MAX_N * MAX_N];	// плоскость столбцевой плотности
	h_print = (box_size/(double)MAX_N);			// шаг вывода

	// цикл формировани€ плоскости
	for(i=0 ; i<MAX_N ; i++)
		for(k=0 ; k<MAX_N ; k++)
		{
			mass_column = 0.0;

			// определение номера €чейки корневой сетки
			i_root = i / MAX_NESTED_N;	
			k_root = k / MAX_NESTED_N;
			
			// суммирование плотности
			for(l_root=0 ; l_root<root_nz ; l_root++)	// цикл по корневой сетке
			{
				l_nested_nz =	// нахождение размера сетки
					Root_Mesh_Hydro[index_root(i_root,k_root,l_root)].nested_N;

				// определение номера €чейки вложенной сетки
				i_nested = (i % MAX_NESTED_N) / (MAX_NESTED_N/l_nested_nz);
				k_nested = (k % MAX_NESTED_N) / (MAX_NESTED_N/l_nested_nz);

				// определение шага вложенной сетки
				h_nested = Root_Mesh_Hydro[index_root(i_root,k_root,l_root)].nested_h;
				
				// цикл по вложенной сетке
				for(l_nested=0 ; l_nested < l_nested_nz ; l_nested++)
				{
					mass_column += Root_Mesh_Hydro[index_root(i_root,k_root,l_root)].
						mesh[Root_Mesh_Hydro[index_root(i_root,k_root,l_root)].
							index_nested_mesh(i_nested,k_nested,l_nested)].rho*
						h_print*h_print*h_nested;
				}
				
			}

			column_density[i*MAX_N+k] = mass_column;
		}

	// генераци€ имени и открытие файла 
	sprintf(FileName,"dens_T%1.3lf.dat",timer);
	fout = fopen(FileName,"w");

	// сохранение плоскости
	total_density = 0.0;
	for(i=0 ; i<MAX_N ; i++)
		for(k=0 ; k<MAX_N ; k++)
		{
			x = i * h_print + h_print/2.0;
			y = k * h_print + h_print/2.0;
			density = column_density[i*MAX_N+k];
			total_density += density;
			fprintf(fout,"%lf %lf %le\n",x,y,density);
		}
	printf("Time %lf\n",timer);

	fclose(fout);
	delete column_density;
}
/*******************************************************************/
#endif