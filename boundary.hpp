#ifndef _boundary_hpp_
#define _boundary_hpp_
/*******************************************************************/
#include "system.hpp"
#include "global.hpp"
#include "riemann.hpp"

// ������� ���������� ��������� ������� �� ���������� ����� �������� �������� �����
void Boundary_Condition()
{
	int i, k, l, im1, km1, lm1, ip1, kp1, lp1, i_astart, k_astart, l_astart;
	int i_nst, k_nst, l_nst, length, i_quad, k_quad, l_quad;
	int nx_nst, ny_nst, nz_nst;
	int nx_nst_m1, nx_nst_p1, ny_nst_m1, ny_nst_p1, nz_nst_m1, nz_nst_p1;
	flux interface_flux, temp_flux, write_flux;
	nested_cell left_cell, right_cell;
	double rho_left, rho_right, velx_left, velx_right, vely_left, vely_right, 
		   velz_left, velz_right, p_left, p_right;
		
	// ���������� ������� ������������� ������� ��� ���� ����� �������� �����
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				// ������ ��������� ����� ������� ������
				nx_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				ny_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				nz_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;

				// ����������� �������� �������� �����
				im1 = (i+root_nx-1) % root_nx;
				km1 = (k+root_ny-1) % root_ny;
				lm1 = (l+root_nz-1) % root_nz;
				ip1 = (i+1) % root_nx;
				kp1 = (k+1) % root_ny;
				lp1 = (l+1) % root_nz;

				// ������ ���������� ��������� ����� �������� �����
				nx_nst_m1 = Root_Mesh_Hydro[index_root(im1,k,l)].nested_N;
				nx_nst_p1 = Root_Mesh_Hydro[index_root(ip1,k,l)].nested_N;
				ny_nst_m1 = Root_Mesh_Hydro[index_root(i,km1,l)].nested_N;
				ny_nst_p1 = Root_Mesh_Hydro[index_root(i,kp1,l)].nested_N;
				nz_nst_m1 = Root_Mesh_Hydro[index_root(i,k,lm1)].nested_N;
				nz_nst_p1 = Root_Mesh_Hydro[index_root(i,k,lp1)].nested_N;

				// ����� ������� ������� �� ��� X
				for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
					for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
					{
						i_nst = 0;				// ������� ���� ����� - ������� 
						i_quad = nx_nst_m1 - 1;	// ��������� ���� ����� - �����

						// ��� ������� ������ ��������� �������
						if( nx_nst < nx_nst_m1 )
						{
							// ��������� ������� �������� ����������
							length = nx_nst_m1/nx_nst;

							// ��������� ��������� ������ ������� �����
							k_astart = k_nst * length;
							l_astart = l_nst * length;

							// ���� ����������� ������� �� ����� ��������
							interface_flux.flux_rho   = 0.0;
							interface_flux.flux_rhovx = 0.0;
							interface_flux.flux_rhovy = 0.0;
							interface_flux.flux_rhovz = 0.0;
							interface_flux.flux_rhoent= 0.0;
							for(k_quad=k_astart ; k_quad<k_astart+length ; k_quad++)
							 for(l_quad=l_astart ; l_quad<l_astart+length ; l_quad++)
							 {
								// ����� ������ 
								left_cell = Root_Mesh_Hydro[index_root(im1,k,l)].
									mesh[Root_Mesh_Hydro[index_root(im1,k,l)].
									index_nested_mesh(i_quad,k_quad,l_quad)]; 
								
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
								temp_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);
								
								// ���������� ������ � ��������� �����
								interface_flux.flux_rho   += temp_flux.flux_rho;
								interface_flux.flux_rhovx += temp_flux.flux_rhovx;
								interface_flux.flux_rhovy += temp_flux.flux_rhovy;
								interface_flux.flux_rhovz += temp_flux.flux_rhovz;
								interface_flux.flux_rhoent+= temp_flux.flux_rhoent;
							 }

							// ������� ������� �� (length*length)
							interface_flux.flux_rho   /= ((double)(length*length));
							interface_flux.flux_rhovx /= ((double)(length*length));
							interface_flux.flux_rhovy /= ((double)(length*length));
							interface_flux.flux_rhovz /= ((double)(length*length));
							interface_flux.flux_rhoent/= ((double)(length*length));
						}
						// ��� ������� ������ ���� �����
						else						
						{
							// ��������� ����� ��������� �� �������
							if(nx_nst == nx_nst_m1)
							{
								k_quad = k_nst;	
								l_quad = l_nst; 
							}

							// �������� ����� ����� ������� ������
							if( nx_nst > nx_nst_m1 )
							{
								length = nx_nst/nx_nst_m1;	// ����������� ������
								k_quad = k_nst/length;		
								l_quad = l_nst/length; 
							}

							// ����� ������ 
							left_cell = Root_Mesh_Hydro[index_root(im1,k,l)].
								mesh[Root_Mesh_Hydro[index_root(im1,k,l)].
								index_nested_mesh(i_quad,k_quad,l_quad)];

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
						}
						// ������ ������� � ������ (ikl_nst) 
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
				   
				// ������ ������� ������� �� ��� X
				for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
					for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
					{
						i_nst = nx_nst-1;	// ��������� ���� ����� - ������� 
						i_quad = 0;			// ������� ���� ����� - ������

						// ��� ������� ������ ��������� �������
						if( nx_nst < nx_nst_p1 )
						{
							// ��������� ������� �������� ����������
							length = nx_nst_p1/nx_nst;

							// ��������� ��������� ������ ������� �����
							k_astart = k_nst * length;
							l_astart = l_nst * length;

							// ���� ����������� ������� �� ����� ��������
							interface_flux.flux_rho   = 0.0;
							interface_flux.flux_rhovx = 0.0;
							interface_flux.flux_rhovy = 0.0;
							interface_flux.flux_rhovz = 0.0;
							interface_flux.flux_rhoent= 0.0;
							for(k_quad=k_astart ; k_quad<k_astart+length ; k_quad++)
							 for(l_quad=l_astart ; l_quad<l_astart+length ; l_quad++)
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
								right_cell = Root_Mesh_Hydro[index_root(ip1,k,l)].
									mesh[Root_Mesh_Hydro[index_root(ip1,k,l)].
									index_nested_mesh(i_quad,k_quad,l_quad)];

								rho_right  = right_cell.rho;
								velx_right = right_cell.velx;
								vely_right = right_cell.vely;
								velz_right = right_cell.velz;
								p_right    = right_cell.press;

								// ������� ������ ������
								temp_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);
								
								// ���������� ������ � ��������� �����
								interface_flux.flux_rho   += temp_flux.flux_rho;
								interface_flux.flux_rhovx += temp_flux.flux_rhovx;
								interface_flux.flux_rhovy += temp_flux.flux_rhovy;
								interface_flux.flux_rhovz += temp_flux.flux_rhovz;
								interface_flux.flux_rhoent+= temp_flux.flux_rhoent;
							 }

							// ������� ������� �� (length*length)
							interface_flux.flux_rho   /= ((double)(length*length));
							interface_flux.flux_rhovx /= ((double)(length*length));
							interface_flux.flux_rhovy /= ((double)(length*length));
							interface_flux.flux_rhovz /= ((double)(length*length));
							interface_flux.flux_rhoent/= ((double)(length*length));
						}
						// ��� ������� ������ ���� �����
						else						
						{
							// ��������� ����� ��������� �� �������
							if(nx_nst == nx_nst_p1)
							{
								k_quad = k_nst;	
								l_quad = l_nst; 
							}

							// �������� ����� ����� ������� ������
							if( nx_nst > nx_nst_p1 )
							{
								length = nx_nst/nx_nst_p1;	// ����������� ������
								k_quad = k_nst/length;		
								l_quad = l_nst/length; 
							}

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
							right_cell = Root_Mesh_Hydro[index_root(ip1,k,l)].
								mesh[Root_Mesh_Hydro[index_root(ip1,k,l)].
								index_nested_mesh(i_quad,k_quad,l_quad)];

							rho_right  = right_cell.rho;
							velx_right = right_cell.velx;
							vely_right = right_cell.vely;
							velz_right = right_cell.velz;
							p_right    = right_cell.press;
							
							// ������� ������ ������
							interface_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);
						}
						// ������ ������� � ������ (ikl_nst) 
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

				// ������ ������� ������� �� ��� Y
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
					{
						k_nst = 0;				// ������� ���� ����� - ������� 
						k_quad = ny_nst_m1 - 1;	// ��������� ���� ����� - �����

						// ��� ������� ������ ��������� �������
						if( ny_nst < ny_nst_m1 )
						{
							// ��������� ������� �������� ����������
							length = ny_nst_m1/ny_nst;

							// ��������� ��������� ������ ������� �����
							i_astart = i_nst * length;
							l_astart = l_nst * length;

							// ���� ����������� ������� �� ����� ��������
							interface_flux.flux_rho   = 0.0;
							interface_flux.flux_rhovx = 0.0;
							interface_flux.flux_rhovy = 0.0;
							interface_flux.flux_rhovz = 0.0;
							interface_flux.flux_rhoent= 0.0;
							for(i_quad=i_astart ; i_quad<i_astart+length ; i_quad++)
							 for(l_quad=l_astart ; l_quad<l_astart+length ; l_quad++)
							 {
								// ����� ������ 
								left_cell = Root_Mesh_Hydro[index_root(i,km1,l)].
									mesh[Root_Mesh_Hydro[index_root(i,km1,l)].
									index_nested_mesh(i_quad,k_quad,l_quad)];

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
								temp_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);
								
								// ���������� ������ � ��������� �����
								interface_flux.flux_rho   += temp_flux.flux_rho;
								interface_flux.flux_rhovx += temp_flux.flux_rhovx;
								interface_flux.flux_rhovy += temp_flux.flux_rhovy;
								interface_flux.flux_rhovz += temp_flux.flux_rhovz;
								interface_flux.flux_rhoent+= temp_flux.flux_rhoent;
							 }

							// ������� ������� �� (length*length)
							interface_flux.flux_rho   /= ((double)(length*length));
							interface_flux.flux_rhovx /= ((double)(length*length));
							interface_flux.flux_rhovy /= ((double)(length*length));
							interface_flux.flux_rhovz /= ((double)(length*length));
							interface_flux.flux_rhoent/= ((double)(length*length));
						}
						// ��� ������� ������ ���� �����
						else						
						{
							// ��������� ����� ��������� �� �������
							if(ny_nst == ny_nst_m1)
							{
								i_quad = i_nst;	
								l_quad = l_nst; 
							}

							// �������� ����� ����� ������� ������
							if( ny_nst > ny_nst_m1 )
							{
								length = ny_nst/ny_nst_m1;	// ����������� ������
								i_quad = i_nst/length;		
								l_quad = l_nst/length; 
							}

							// ����� ������ 
							left_cell = Root_Mesh_Hydro[index_root(i,km1,l)].
								mesh[Root_Mesh_Hydro[index_root(i,km1,l)].
								index_nested_mesh(i_quad,k_quad,l_quad)];

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
						}
						// ������ ������� � ������ (ikl_nst) 
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

				// ������� ������� ������� �� ��� Y
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(l_nst=0 ; l_nst<nz_nst ; l_nst++)
					{
						k_nst = ny_nst-1;	// ��������� ���� ����� - ������� 
						k_quad = 0;			// ������� ���� ����� - ������

						// ��� ������� ������ ��������� �������
						if( ny_nst < ny_nst_p1 )
						{
							// ��������� ������� �������� ����������
							length = ny_nst_p1/ny_nst;

							// ��������� ��������� ������ ������� �����
							i_astart = i_nst * length;
							l_astart = l_nst * length;

							// ���� ����������� ������� �� ����� ��������
							interface_flux.flux_rho   = 0.0;
							interface_flux.flux_rhovx = 0.0;
							interface_flux.flux_rhovy = 0.0;
							interface_flux.flux_rhovz = 0.0;
							interface_flux.flux_rhoent= 0.0;
							for(i_quad=i_astart ; i_quad<i_astart+length ; i_quad++)
							 for(l_quad=l_astart ; l_quad<l_astart+length ; l_quad++)
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
								right_cell = Root_Mesh_Hydro[index_root(i,kp1,l)].
									mesh[Root_Mesh_Hydro[index_root(i,kp1,l)].
									index_nested_mesh(i_quad,k_quad,l_quad)];

								rho_right  = right_cell.rho;
								velx_right = right_cell.vely;
								vely_right = right_cell.velz;
								velz_right = right_cell.velx;
								p_right    = right_cell.press;

								// ������� ������ ������
								temp_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);
								
								// ���������� ������ � ��������� �����
								interface_flux.flux_rho   += temp_flux.flux_rho;
								interface_flux.flux_rhovx += temp_flux.flux_rhovx;
								interface_flux.flux_rhovy += temp_flux.flux_rhovy;
								interface_flux.flux_rhovz += temp_flux.flux_rhovz;
								interface_flux.flux_rhoent+= temp_flux.flux_rhoent;
							 }

							// ������� ������� �� (length*length)
							interface_flux.flux_rho   /= ((double)(length*length));
							interface_flux.flux_rhovx /= ((double)(length*length));
							interface_flux.flux_rhovy /= ((double)(length*length));
							interface_flux.flux_rhovz /= ((double)(length*length));
							interface_flux.flux_rhoent/= ((double)(length*length));
						}
						// ��� ������� ������ ���� �����
						else						
						{
							// ��������� ����� ��������� �� �������
							if(ny_nst == ny_nst_p1)
							{
								i_quad = i_nst;	
								l_quad = l_nst; 
							}

							// �������� ����� ����� ������� ������
							if( ny_nst > ny_nst_p1 )
							{
								length = ny_nst/ny_nst_p1;	// ����������� ������
								i_quad = i_nst/length;		
								l_quad = l_nst/length; 
							}

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
							right_cell = Root_Mesh_Hydro[index_root(i,kp1,l)].
								mesh[Root_Mesh_Hydro[index_root(i,kp1,l)].
								index_nested_mesh(i_quad,k_quad,l_quad)];

							rho_right  = right_cell.rho;
							velx_right = right_cell.vely;
							vely_right = right_cell.velz;
							velz_right = right_cell.velx;
							p_right    = right_cell.press;
							
							// ������� ������ ������
							interface_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);
						}
						// ������ ������� � ������ (ikl_nst) 
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

				// ������� ������� ������� �� ��� Z
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
					{
						l_nst = 0;				// ������� ���� ����� - ������� 
						l_quad = nz_nst_m1 - 1;	// ��������� ���� ����� - �����

						// ��� ������� ������ ��������� �������
						if( nz_nst < nz_nst_m1 )
						{
							// ��������� ������� �������� ����������
							length = nz_nst_m1/nz_nst;

							// ��������� ��������� ������ ������� �����
							i_astart = i_nst * length;
							k_astart = k_nst * length;

							// ���� ����������� ������� �� ����� ��������
							interface_flux.flux_rho   = 0.0;
							interface_flux.flux_rhovx = 0.0;
							interface_flux.flux_rhovy = 0.0;
							interface_flux.flux_rhovz = 0.0;
							interface_flux.flux_rhoent= 0.0;
							for(i_quad=i_astart ; i_quad<i_astart+length ; i_quad++)
							 for(k_quad=k_astart ; k_quad<k_astart+length ; k_quad++)
							 {
								// ����� ������ 
							    left_cell = Root_Mesh_Hydro[index_root(i,k,lm1)].
									mesh[Root_Mesh_Hydro[index_root(i,k,lm1)].
									index_nested_mesh(i_quad,k_quad,l_quad)];

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
								temp_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);
								
								// ���������� ������ � ��������� �����
								interface_flux.flux_rho   += temp_flux.flux_rho;
								interface_flux.flux_rhovx += temp_flux.flux_rhovx;
								interface_flux.flux_rhovy += temp_flux.flux_rhovy;
								interface_flux.flux_rhovz += temp_flux.flux_rhovz;
								interface_flux.flux_rhoent+= temp_flux.flux_rhoent;
							 }

							// ������� ������� �� (length*length)
							interface_flux.flux_rho   /= ((double)(length*length));
							interface_flux.flux_rhovx /= ((double)(length*length));
							interface_flux.flux_rhovy /= ((double)(length*length));
							interface_flux.flux_rhovz /= ((double)(length*length));
							interface_flux.flux_rhoent/= ((double)(length*length));
						}
						// ��� ������� ������ ���� �����
						else						
						{
							// ��������� ����� ��������� �� �������
							if(nz_nst == nz_nst_m1)
							{
								i_quad = i_nst;	
								k_quad = k_nst; 
							}

							// �������� ����� ����� ������� ������
							if( nz_nst > nz_nst_m1 )
							{
								length = nz_nst/nz_nst_m1;	// ����������� ������
								i_quad = i_nst/length;		
								k_quad = k_nst/length; 
							}

							// ����� ������ 
							left_cell = Root_Mesh_Hydro[index_root(i,k,lm1)].
								mesh[Root_Mesh_Hydro[index_root(i,k,lm1)].
								index_nested_mesh(i_quad,k_quad,l_quad)];

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
						}
						// ������ ������� � ������ (ikl_nst) 
						write_flux.flux_rho = interface_flux.flux_rho;
						write_flux.flux_rhovx = interface_flux.flux_rhovy;
						write_flux.flux_rhovy = interface_flux.flux_rhovz;
						write_flux.flux_rhovz = interface_flux.flux_rhovx;
						write_flux.flux_rhoent= interface_flux.flux_rhoent;

						Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].
								z_minus = write_flux;
					}

				// ������� ������� ������� �� ��� Z
				for(i_nst=0 ; i_nst<nx_nst ; i_nst++)
					for(k_nst=0 ; k_nst<ny_nst ; k_nst++)
					{
						l_nst = nz_nst-1;	// ��������� ���� ����� - ������� 
						l_quad = 0;			// ������� ���� ����� - ������

						// ��� ������� ������ ��������� �������
						if( nz_nst < nz_nst_p1 )
						{
							// ��������� ������� �������� ����������
							length = nz_nst_p1/nz_nst;

							// ��������� ��������� ������ ������� �����
							i_astart = i_nst * length;
							k_astart = k_nst * length;

							// ���� ����������� ������� �� ����� ��������
							interface_flux.flux_rho   = 0.0;
							interface_flux.flux_rhovx = 0.0;
							interface_flux.flux_rhovy = 0.0;
							interface_flux.flux_rhovz = 0.0;
							interface_flux.flux_rhoent= 0.0;
							for(i_quad=i_astart ; i_quad<i_astart+length ; i_quad++)
							 for(k_quad=k_astart ; k_quad<k_astart+length ; k_quad++)
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
								right_cell = Root_Mesh_Hydro[index_root(i,k,lp1)].
									mesh[Root_Mesh_Hydro[index_root(i,k,lp1)].
									index_nested_mesh(i_quad,k_quad,l_quad)];

								rho_right  = right_cell.rho;
								velx_right = right_cell.velz;
								vely_right = right_cell.velx;
								velz_right = right_cell.vely;
								p_right    = right_cell.press;

								// ������� ������ ������
								temp_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);
								
								// ���������� ������ � ��������� �����
								interface_flux.flux_rho   += temp_flux.flux_rho;
								interface_flux.flux_rhovx += temp_flux.flux_rhovx;
								interface_flux.flux_rhovy += temp_flux.flux_rhovy;
								interface_flux.flux_rhovz += temp_flux.flux_rhovz;
								interface_flux.flux_rhoent+= temp_flux.flux_rhoent;
							 }

							// ������� ������� �� (length*length)
							interface_flux.flux_rho   /= ((double)(length*length));
							interface_flux.flux_rhovx /= ((double)(length*length));
							interface_flux.flux_rhovy /= ((double)(length*length));
							interface_flux.flux_rhovz /= ((double)(length*length));
							interface_flux.flux_rhoent/= ((double)(length*length));
						}
						// ��� ������� ������ ���� �����
						else						
						{
							// ��������� ����� ��������� �� �������
							if(nz_nst == nz_nst_p1)
							{
								i_quad = i_nst;	
								k_quad = k_nst; 
							}

							// �������� ����� ����� ������� ������
							if( nz_nst > nz_nst_p1 )
							{
								length = nz_nst/nz_nst_p1;	// ����������� ������
								i_quad = i_nst/length;		
								k_quad = k_nst/length; 
							}

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
							right_cell = Root_Mesh_Hydro[index_root(i,k,lp1)].
								mesh[Root_Mesh_Hydro[index_root(i,k,lp1)].
								index_nested_mesh(i_quad,k_quad,l_quad)];

							rho_right  = right_cell.rho;
							velx_right = right_cell.velz;
							vely_right = right_cell.velx;
							velz_right = right_cell.vely;
							p_right    = right_cell.press;
							
							// ������� ������ ������
							interface_flux = riemann(rho_left, rho_right, 
									velx_left, velx_right, vely_left, vely_right, 
									velz_left, velz_right, p_left, p_right);
						}
						// ������ ������� � ������ (ikl_nst) 
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
/*******************************************************************/
#endif