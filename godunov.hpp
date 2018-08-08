#ifndef _godunov_hpp_
#define _godunov_hpp_
/*******************************************************************/
#include "system.hpp"
#include "global.hpp"
#include "riemann.hpp"

// ������� ���������� ������� ��� ������ ��������� �����
void Riemann_Computing()
{
	int i, k, l, nx_nst, ny_nst, nz_nst, i_nst, k_nst, l_nst;
	flux interface_flux, write_flux;
	nested_cell left_cell, right_cell;
	double rho_left, rho_right, velx_left, velx_right, vely_left, vely_right, 
		   velz_left, velz_right, p_left, p_right;

	// ���������� ������� ������ ���� ��������� �����
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				nx_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				ny_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				nz_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;

				// ������ � ��������� ������
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
						for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
						{
							// ����� ��������� �� ��� X
							if(i_nst > 0)
							{
								// ����� ������ 
								left_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst-1,k_nst,l_nst)];

								rho_left  = left_cell.rho;
								velx_left = left_cell.velx;
								vely_left = left_cell.vely;
								velz_left = left_cell.velz;
								p_left    = left_cell.press;

								// ������ ������
								right_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)];

								rho_right  = right_cell.rho;
								velx_right = right_cell.velx;
								vely_right = right_cell.vely;
								velz_right = right_cell.velz;
								p_right    = right_cell.press;

								// ������� ������ ������
								interface_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);

								// ������ ������� � ������
								write_flux.flux_rho   = interface_flux.flux_rho;
								write_flux.flux_rhovx = interface_flux.flux_rhovx;
								write_flux.flux_rhovy = interface_flux.flux_rhovy;
								write_flux.flux_rhovz = interface_flux.flux_rhovz;
								write_flux.flux_rhoent= interface_flux.flux_rhoent;

								Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)].
									x_minus = write_flux;
							}

							// ������ ��������� �� ��� X
							if(i_nst < nx_nst-1)
							{
								// ����� ������ 
								left_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)];

								rho_left  = left_cell.rho;
								velx_left = left_cell.velx;
								vely_left = left_cell.vely;
								velz_left = left_cell.velz;
								p_left    = left_cell.press;

								// ������ ������
								right_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst+1,k_nst,l_nst)];

								rho_right  = right_cell.rho;
								velx_right = right_cell.velx;
								vely_right = right_cell.vely;
								velz_right = right_cell.velz;
								p_right    = right_cell.press;

								// ������� ������ ������
								interface_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);

								// ������ ������� � ������
								write_flux.flux_rho   = interface_flux.flux_rho;
								write_flux.flux_rhovx = interface_flux.flux_rhovx;
								write_flux.flux_rhovy = interface_flux.flux_rhovy;
								write_flux.flux_rhovz = interface_flux.flux_rhovz;
								write_flux.flux_rhoent= interface_flux.flux_rhoent;

								Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)].
									x_plus = write_flux;
							}

							// ������ ��������� �� ��� Y
							if(k_nst > 0)
							{
								// ����� ������ 
								left_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst-1,l_nst)];

								rho_left  = left_cell.rho;
								velx_left = left_cell.vely;
								vely_left = left_cell.velz;
								velz_left = left_cell.velx;
								p_left    = left_cell.press;

								// ������ ������
								right_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)];

								rho_right  = right_cell.rho;
								velx_right = right_cell.vely;
								vely_right = right_cell.velz;
								velz_right = right_cell.velx;
								p_right    = right_cell.press;

								// ������� ������ ������
								interface_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);

								// ������ ������� � ������
								write_flux.flux_rho   = interface_flux.flux_rho;
								write_flux.flux_rhovx = interface_flux.flux_rhovz;
								write_flux.flux_rhovy = interface_flux.flux_rhovx;
								write_flux.flux_rhovz = interface_flux.flux_rhovy;
								write_flux.flux_rhoent= interface_flux.flux_rhoent;

								Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)].
									y_minus = write_flux;
							}

							// ������� ��������� �� ��� Y
							if(k_nst < ny_nst-1)
							{
								// ����� ������ 
								left_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)];

								rho_left  = left_cell.rho;
								velx_left = left_cell.vely;
								vely_left = left_cell.velz;
								velz_left = left_cell.velx;
								p_left    = left_cell.press;

								// ������ ������
								right_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst+1,l_nst)];

								rho_right  = right_cell.rho;
								velx_right = right_cell.vely;
								vely_right = right_cell.velz;
								velz_right = right_cell.velx;
								p_right    = right_cell.press;

								// ������� ������ ������
								interface_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);

								// ������ ������� � ������
								write_flux.flux_rho   = interface_flux.flux_rho;
								write_flux.flux_rhovx = interface_flux.flux_rhovz;
								write_flux.flux_rhovy = interface_flux.flux_rhovx;
								write_flux.flux_rhovz = interface_flux.flux_rhovy;
								write_flux.flux_rhoent= interface_flux.flux_rhoent;

								Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)].
									y_plus = write_flux;
							}

							// ������� ��������� �� ��� Z
							if(l_nst > 0)
							{
								// ����� ������ 
								left_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst-1)];

								rho_left  = left_cell.rho;
								velx_left = left_cell.velz;
								vely_left = left_cell.velx;
								velz_left = left_cell.vely;
								p_left    = left_cell.press;

								// ������ ������
								right_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)];

								rho_right  = right_cell.rho;
								velx_right = right_cell.velz;
								vely_right = right_cell.velx;
								velz_right = right_cell.vely;
								p_right    = right_cell.press;

								// ������� ������ ������
								interface_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);

								// ������ ������� � ������
								write_flux.flux_rho   = interface_flux.flux_rho;
								write_flux.flux_rhovx = interface_flux.flux_rhovy;
								write_flux.flux_rhovy = interface_flux.flux_rhovz;
								write_flux.flux_rhovz = interface_flux.flux_rhovx;
								write_flux.flux_rhoent= interface_flux.flux_rhoent;

								Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)].
									z_minus = write_flux;
							}

							// ������� �������� �� ��� Z
							if(l_nst < nz_nst-1)
							{
								// ����� ������ 
								left_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)];

								rho_left  = left_cell.rho;
								velx_left = left_cell.velz;
								vely_left = left_cell.velx;
								velz_left = left_cell.vely;
								p_left    = left_cell.press;

								// ������ ������
								right_cell = Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst+1)];

								rho_right  = right_cell.rho;
								velx_right = right_cell.velz;
								vely_right = right_cell.velx;
								velz_right = right_cell.vely;
								p_right    = right_cell.press;

								// ������� ������ ������
								interface_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);

								// ������ ������� � ������
								write_flux.flux_rho   = interface_flux.flux_rho;
								write_flux.flux_rhovx = interface_flux.flux_rhovy;
								write_flux.flux_rhovy = interface_flux.flux_rhovz;
								write_flux.flux_rhovz = interface_flux.flux_rhovx;
								write_flux.flux_rhoent= interface_flux.flux_rhoent;

								Root_Mesh_Hydro[index_root(i,k,l)].
									mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_mesh(i_nst,k_nst,l_nst)].
									z_plus = write_flux;
							}
						}
			}
}

// ������� ������� ��������� ������������� ������� ��������
void Godunov_Solver()
{
	int i, k, l, nx_nst, ny_nst, nz_nst, i_nst, k_nst, l_nst;
	double h_nst;
	nested_cell cell;

	// ���������� ����� ��������
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
							// ������ ������
							cell = Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)];

							// ���������� ����� �������� � ������
							// ���������
							cell.rho = cell.rho - tau * (
								(cell.x_plus.flux_rho-cell.x_minus.flux_rho)/h_nst + 
								(cell.y_plus.flux_rho-cell.y_minus.flux_rho)/h_nst +
								(cell.z_plus.flux_rho-cell.z_minus.flux_rho)/h_nst);

							// x-���������� ��������
							cell.rhovx = cell.rhovx - tau * (
								(cell.x_plus.flux_rhovx-cell.x_minus.flux_rhovx)/h_nst + 
								(cell.y_plus.flux_rhovx-cell.y_minus.flux_rhovx)/h_nst +
								(cell.z_plus.flux_rhovx-cell.z_minus.flux_rhovx)/h_nst);

							// y-���������� ��������
							cell.rhovy = cell.rhovy - tau * (
								(cell.x_plus.flux_rhovy-cell.x_minus.flux_rhovy)/h_nst + 
								(cell.y_plus.flux_rhovy-cell.y_minus.flux_rhovy)/h_nst +
								(cell.z_plus.flux_rhovy-cell.z_minus.flux_rhovy)/h_nst);

							// z-���������� ��������
							cell.rhovz = cell.rhovz - tau * (
								(cell.x_plus.flux_rhovz-cell.x_minus.flux_rhovz)/h_nst + 
								(cell.y_plus.flux_rhovz-cell.y_minus.flux_rhovz)/h_nst +
								(cell.z_plus.flux_rhovz-cell.z_minus.flux_rhovz)/h_nst);

							// ��������� ��������
							cell.rhoent = cell.rhoent - tau * (
								(cell.x_plus.flux_rhoent-cell.x_minus.flux_rhoent)/h_nst + 
								(cell.y_plus.flux_rhoent-cell.y_minus.flux_rhoent)/h_nst +
								(cell.z_plus.flux_rhoent-cell.z_minus.flux_rhoent)/h_nst);

							// ������ ������
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)] = cell;
						}
			}
}
/*******************************************************************/
#endif