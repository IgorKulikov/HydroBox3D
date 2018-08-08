#ifndef _setup_hpp_
#define _setup_hpp_
/*******************************************************************/
#include "global.hpp"
#include "system.hpp"
#include "physics.hpp"

/* ������� �������� �������� � ��������� */
void Init()
{
	Root_Mesh_Hydro = new root_cell[root_nx*root_ny*root_nz];
	xfft = new double[root_nx];
	yfft = new double[root_nx];
	fout = fopen("output.log","w");
}

/* ������� ������������ ������ */
void Setup()
{
	int i, k, l, i_nst, k_nst, l_nst, nx_nst, ny_nst, nz_nst;
	double x, y, z, x0, y0, z0, rad, h_nst;
	double density_setup, pressure_setup, entropy_setup;
	
	// ����������� ������ �������
	x0 = box_size/2.0;
	y0 = box_size/2.0;
	z0 = box_size/2.0;

	// �������� ��������� ��������� �����
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				x = i*root_h + root_h/2.0 - x0;
				y = k*root_h + root_h/2.0 - y0;
				z = l*root_h + root_h/2.0 - z0;
				rad = sqrt(x*x+y*y+z*z);

				// [8x8x8] 0.5 [4x4x4] 1.0 [2x2x2]
				if(rad <= 0.5)	
					Root_Mesh_Hydro[index_root(i,k,l)].make_nested_mesh(8);
				else
					if(rad <= 1.0)
						Root_Mesh_Hydro[index_root(i,k,l)].make_nested_mesh(4);
					else
						Root_Mesh_Hydro[index_root(i,k,l)].make_nested_mesh(2);

			}

	// ������� ��������� ������
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				nx_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				ny_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				nz_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				h_nst  = Root_Mesh_Hydro[index_root(i,k,l)].nested_h;

				// ������ � ��������� ������
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
						for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
						{
							
							x = i*root_h + i_nst*h_nst + h_nst/2.0 - x0;
							y = k*root_h + k_nst*h_nst + h_nst/2.0 - y0;
							z = l*root_h + l_nst*h_nst + h_nst/2.0 - z0;
							rad = sqrt(x*x+y*y+z*z);

							// ���������� ��������� � ��������
							density_setup = density_profile(rad); 
							pressure_setup = pressure_profile(rad);
							entropy_setup = get_entropy(density_setup,pressure_setup);

							// ���������
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rho = 
									density_setup;

							// ������ ��������
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhovx = 0.0;
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhovy = 0.0;
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhovz = 0.0;

							// ��������� ��������
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhoent = 
									density_setup * entropy_setup;
						}
			}
}

/* ������� �������� �������� � ��������� */
void Destroy()
{
	int i, k, l;

	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
				Root_Mesh_Hydro[index_root(i,k,l)].destroy_nested_mesh();
			
	delete Root_Mesh_Hydro;

	delete xfft;
	delete yfft;

	fclose(fout);
}
/*******************************************************************/
#endif