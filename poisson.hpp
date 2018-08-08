#ifndef _poisson_hpp_
#define _poisson_hpp_
/*******************************************************************/
#include <math.h>
#include <memory.h>
#include "system.hpp"
#include "fft.hpp"

// �������� ���������� �� �������� ������ � ����
double gravity_cell2node(int i, int ipm, int k, int kpm, int l, int lpm)
{	
	return 0.125*(	Root_Mesh_Hydro[index_root(i,  k,  l  )].phi + 
					Root_Mesh_Hydro[index_root(i,  k,  lpm)].phi + 
					Root_Mesh_Hydro[index_root(i,  kpm,l  )].phi + 
					Root_Mesh_Hydro[index_root(i,  kpm,lpm)].phi + 
					Root_Mesh_Hydro[index_root(ipm,k,  l  )].phi + 
					Root_Mesh_Hydro[index_root(ipm,k,  lpm)].phi + 
					Root_Mesh_Hydro[index_root(ipm,kpm,l  )].phi + 
					Root_Mesh_Hydro[index_root(ipm,kpm,lpm)].phi );
}

// ������� �������� �������������� �����
void FFT3D(int dir)
{
	int i, k, l;

	// �������������� �� ��� X
	for(k=0 ; k<root_ny ; k++)
		for(l=0 ; l<root_nz ; l++)
			{
				for(i=0 ; i<root_nx ; i++)
				{
					xfft[i] = Root_Mesh_Hydro[index_root(i,k,l)].phi_real;
					yfft[i] = Root_Mesh_Hydro[index_root(i,k,l)].phi_image;
				}
				
				FFT(dir,PSN_POW,xfft,yfft);
				
				for(i=0 ; i<root_nx ; i++)
				{
					Root_Mesh_Hydro[index_root(i,k,l)].phi_real  = xfft[i];
					Root_Mesh_Hydro[index_root(i,k,l)].phi_image = yfft[i];
				}
			}
	
	// �������������� �� ��� Y
	for(i=0 ; i<root_nx ; i++)
		for(l=0 ; l<root_nz ; l++)
			{
				for(k=0 ; k<root_ny ; k++)
				{
					xfft[k] = Root_Mesh_Hydro[index_root(i,k,l)].phi_real;
					yfft[k] = Root_Mesh_Hydro[index_root(i,k,l)].phi_image;
				}
				
				FFT(dir,PSN_POW,xfft,yfft);
				
				for(k=0 ; k<root_ny ; k++)
				{
					Root_Mesh_Hydro[index_root(i,k,l)].phi_real  = xfft[k];
					Root_Mesh_Hydro[index_root(i,k,l)].phi_image = yfft[k];
				}
			}

	// �������������� �� ��� Z
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			{
				for(l=0 ; l<root_nz ; l++)
				{
					xfft[l] = Root_Mesh_Hydro[index_root(i,k,l)].phi_real;
					yfft[l] = Root_Mesh_Hydro[index_root(i,k,l)].phi_image;
				}
				
				FFT(dir,PSN_POW,xfft,yfft);
				
				for(l=0 ; l<root_nz ; l++)
				{
					Root_Mesh_Hydro[index_root(i,k,l)].phi_real  = xfft[l];
					Root_Mesh_Hydro[index_root(i,k,l)].phi_image = yfft[l];
				}
			}
}

// ������� ��������� ��������
void Poisson_Solver()
{
	int i, k, l, i_nst, k_nst, l_nst, nx_nst, ny_nst, nz_nst;
	int im1, ip1, km1, kp1, lp1, lm1, sor_num_iter;
	double phi_ip_kp_lp, phi_ip_kp_lm, phi_ip_km_lp, phi_ip_km_lm,
		   phi_im_kp_lp, phi_im_kp_lm, phi_im_km_lp, phi_im_km_lm;
	double x, y, z, x0, y0, z0, rad, h_nst, core_mass, residual, med_density;
	double nested_phi, coeff, schemas_i, schemas_k, schemas_l;

	// ������� ��������� �������� �� �������� �����
	core_mass = 0.0;
	
	// ������� ������ �����
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				core_mass += 
					Root_Mesh_Hydro[index_root(i,k,l)].rho * root_h * root_h * root_h;

				Root_Mesh_Hydro[index_root(i,k,l)].phi_real = 
					4.0*PI*Root_Mesh_Hydro[index_root(i,k,l)].rho;

				Root_Mesh_Hydro[index_root(i,k,l)].phi_image = 0.0;
			}
	
	// ������ �������������� �����
	FFT3D(1);

	// ���������� �����
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				schemas_i = 1.0-(2.0/3.0)*(sin((i*PI)/root_nx)*sin((i*PI)/root_nx));
				schemas_k = 1.0-(2.0/3.0)*(sin((k*PI)/root_nx)*sin((k*PI)/root_nx));
				schemas_l = 1.0-(2.0/3.0)*(sin((l*PI)/root_nx)*sin((l*PI)/root_nx));

				if( i+k+l == 0 ) 
					coeff = 0.0;
				else 
					coeff = -root_h*root_h/(6.0*(1.0-schemas_i*schemas_k*schemas_l));

				Root_Mesh_Hydro[index_root(i,k,l)].phi_real  *= coeff;
				Root_Mesh_Hydro[index_root(i,k,l)].phi_image *= coeff;
			}
	
	// �������� �������������� �����
	FFT3D(-1);

	// ������ ����������
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
				Root_Mesh_Hydro[index_root(i,k,l)].phi = 
					Root_Mesh_Hydro[index_root(i,k,l)].phi_real;
	
	/* ������� ���������� ���������� �������� ����� 
	x0 = box_size/2.0;		// ����������� ������ �������
	y0 = box_size/2.0;
	z0 = box_size/2.0;
	core_mass = 1.0;		// �������� �����

	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				x = i*root_h + root_h/2.0 - x0;
				y = k*root_h + root_h/2.0 - y0;
				z = l*root_h + root_h/2.0 - z0;
				rad = sqrt(x*x+y*y+z*z);
				Root_Mesh_Hydro[index_root(i,k,l)].phi = -core_mass/rad;
			}
	*/

	// �������� ���������� �������� ����� �� ��������� �����
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				// ������ ���������� ��������� �����
				nx_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				ny_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				nz_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				h_nst  = Root_Mesh_Hydro[index_root(i,k,l)].nested_h;

				// ����������� �������� �������� �����
				im1 = (i+root_nx-1) % root_nx;
				km1 = (k+root_ny-1) % root_ny;
				lm1 = (l+root_nz-1) % root_nz;
				ip1 = (i+1) % root_nx;
				kp1 = (k+1) % root_ny;
				lp1 = (l+1) % root_nz;

				// ����������� ������� ����� ����������
				phi_ip_kp_lp = gravity_cell2node(i,ip1,k,kp1,l,lp1);
				phi_ip_kp_lm = gravity_cell2node(i,ip1,k,kp1,l,lm1);
				phi_ip_km_lp = gravity_cell2node(i,ip1,k,km1,l,lp1);
				phi_ip_km_lm = gravity_cell2node(i,ip1,k,km1,l,lm1);
				phi_im_kp_lp = gravity_cell2node(i,im1,k,kp1,l,lp1);
				phi_im_kp_lm = gravity_cell2node(i,im1,k,kp1,l,lm1);
				phi_im_km_lp = gravity_cell2node(i,im1,k,km1,l,lp1);
				phi_im_km_lm = gravity_cell2node(i,im1,k,km1,l,lm1);

				// �������� ���������� �� ��������� �����
				for(i_nst=0 ; i_nst<=nx_nst ; i_nst++)
					for(k_nst=0 ; k_nst<=ny_nst ; k_nst++)
						for(l_nst=0 ; l_nst<=nz_nst ; l_nst++)
						{
							/* ������� �������� ����������� �� ����� */
							x = i_nst * h_nst;	// ���������� ��������� �����
							y = k_nst * h_nst;
							z = l_nst * h_nst;

							nested_phi =	// ���������� ����� �� �����������
							 phi_ip_kp_lp*(0.0+x/root_h)*(0.0+y/root_h)*(0.0+z/root_h)+ 
							 phi_ip_kp_lm*(0.0+x/root_h)*(0.0+y/root_h)*(1.0-z/root_h)+ 
							 phi_ip_km_lp*(0.0+x/root_h)*(1.0-y/root_h)*(0.0+z/root_h)+ 
							 phi_ip_km_lm*(0.0+x/root_h)*(1.0-y/root_h)*(1.0-z/root_h)+ 
							 phi_im_kp_lp*(1.0-x/root_h)*(0.0+y/root_h)*(0.0+z/root_h)+ 
							 phi_im_kp_lm*(1.0-x/root_h)*(0.0+y/root_h)*(1.0-z/root_h)+ 
							 phi_im_km_lp*(1.0-x/root_h)*(1.0-y/root_h)*(0.0+z/root_h)+ 
							 phi_im_km_lm*(1.0-x/root_h)*(1.0-y/root_h)*(1.0-z/root_h);

							Root_Mesh_Hydro[index_root(i,k,l)].
								gravity[Root_Mesh_Hydro[index_root(i,k,l)].
									index_nested_gravity_mesh(i_nst,k_nst,l_nst)] = 
										nested_phi;
								
						}
			}

	// ������� ��������� �������� �� ��������� ������
	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				// ������ ���������� ��������� �����
				nx_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				ny_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				nz_nst = Root_Mesh_Hydro[index_root(i,k,l)].nested_N;
				h_nst  = Root_Mesh_Hydro[index_root(i,k,l)].nested_h;

				// ������������ �������
				sor_num_iter = 0;

				// ������������� ��������� �������
				memcpy(	Root_Mesh_Hydro[index_root(i,k,l)].grav_new,
						Root_Mesh_Hydro[index_root(i,k,l)].gravity,
						(nx_nst+1)*(ny_nst+1)*(nz_nst+1)*sizeof(double));

				while(sor_num_iter < MAX_SOR_ITER)
				{
					// ���� �� ���������� ����� ��������� �����
					for(i_nst=1 ; i_nst<nx_nst ; i_nst++)
						for(k_nst=1 ; k_nst<ny_nst ; k_nst++)
							for(l_nst=1 ; l_nst<nz_nst ; l_nst++)
							{
								// ��������� ��������� � ����
								med_density = 0.125 * (
									Root_Mesh_Hydro[index_root(i,k,l)].
									 mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									 index_nested_mesh(i_nst  ,k_nst  ,l_nst  )].rho + 
									Root_Mesh_Hydro[index_root(i,k,l)].
									 mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									 index_nested_mesh(i_nst  ,k_nst  ,l_nst-1)].rho + 
									Root_Mesh_Hydro[index_root(i,k,l)].
									 mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									 index_nested_mesh(i_nst  ,k_nst-1,l_nst  )].rho + 
									Root_Mesh_Hydro[index_root(i,k,l)].
									 mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									 index_nested_mesh(i_nst  ,k_nst-1,l_nst-1)].rho + 
									Root_Mesh_Hydro[index_root(i,k,l)].
									 mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									 index_nested_mesh(i_nst-1,k_nst  ,l_nst  )].rho + 
									Root_Mesh_Hydro[index_root(i,k,l)].
									 mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									 index_nested_mesh(i_nst-1,k_nst  ,l_nst-1)].rho + 
									Root_Mesh_Hydro[index_root(i,k,l)].
									 mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									 index_nested_mesh(i_nst-1,k_nst-1,l_nst  )].rho + 
									Root_Mesh_Hydro[index_root(i,k,l)].
									 mesh[Root_Mesh_Hydro[index_root(i,k,l)].
									 index_nested_mesh(i_nst-1,k_nst-1,l_nst-1)].rho);

								// ��������� �������
								residual = (
								 Root_Mesh_Hydro[index_root(i,k,l)].
								  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
								  index_nested_gravity_mesh(i_nst+1,k_nst  ,l_nst  )]+ 
								 Root_Mesh_Hydro[index_root(i,k,l)].
								  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
								  index_nested_gravity_mesh(i_nst-1,k_nst  ,l_nst  )]+ 
								 Root_Mesh_Hydro[index_root(i,k,l)].
								  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
								  index_nested_gravity_mesh(i_nst  ,k_nst+1,l_nst  )]+ 
								 Root_Mesh_Hydro[index_root(i,k,l)].
								  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
								  index_nested_gravity_mesh(i_nst  ,k_nst-1,l_nst  )]+ 
								 Root_Mesh_Hydro[index_root(i,k,l)].
								  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
								  index_nested_gravity_mesh(i_nst  ,k_nst  ,l_nst+1)]+ 
								 Root_Mesh_Hydro[index_root(i,k,l)].
								  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
								  index_nested_gravity_mesh(i_nst  ,k_nst  ,l_nst-1)]-
								 6.0*Root_Mesh_Hydro[index_root(i,k,l)].
								  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
								  index_nested_gravity_mesh(i_nst,k_nst,l_nst)]) - 
								  4.0*PI*med_density*h_nst*h_nst;

								// ��������� ����� ����������� ����������
								Root_Mesh_Hydro[index_root(i,k,l)].
								  grav_new[Root_Mesh_Hydro[index_root(i,k,l)].
								  index_nested_gravity_mesh(i_nst,k_nst,l_nst)] =
								Root_Mesh_Hydro[index_root(i,k,l)].
								  gravity[Root_Mesh_Hydro[index_root(i,k,l)].
								  index_nested_gravity_mesh(i_nst,k_nst,l_nst)] +
								RELAXATION_SOR*residual/6.0;
							}

					// ��������� ��������� ����������� ����������
					memcpy(	Root_Mesh_Hydro[index_root(i,k,l)].gravity,
							Root_Mesh_Hydro[index_root(i,k,l)].grav_new,
							(nx_nst+1)*(ny_nst+1)*(nz_nst+1)*sizeof(double));

					// ��������� � ��������� ��������
					sor_num_iter++;
				}
			}
}

// ���������� �������������� ����
void Gravity_Force()
{	
	int i, k, l, i_nst, k_nst, l_nst, nx_nst, ny_nst, nz_nst;
	double h_nst, gradfix, gradfiy, gradfiz;
	double phi_ip_kp_lp, phi_ip_kp_lm, phi_ip_km_lp, phi_ip_km_lm,
		   phi_im_kp_lp, phi_im_kp_lm, phi_im_km_lp, phi_im_km_lm;
	nested_cell cell;

	for(i=0 ; i<root_nx ; i++)
		for(k=0 ; k<root_ny ; k++)
			for(l=0 ; l<root_nz ; l++)
			{
				// ������ ���������� �������� ������
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

							// ������ ���������� � �����
							phi_ip_kp_lp = Root_Mesh_Hydro[index_root(i,k,l)].gravity
								[Root_Mesh_Hydro[index_root(i,k,l)].
								 index_nested_gravity_mesh(i_nst+1,k_nst+1,l_nst+1)];
							phi_ip_kp_lm = Root_Mesh_Hydro[index_root(i,k,l)].gravity
								[Root_Mesh_Hydro[index_root(i,k,l)].
								 index_nested_gravity_mesh(i_nst+1,k_nst+1,l_nst  )];
							phi_ip_km_lp = Root_Mesh_Hydro[index_root(i,k,l)].gravity
								[Root_Mesh_Hydro[index_root(i,k,l)].
								 index_nested_gravity_mesh(i_nst+1,k_nst  ,l_nst+1)];
							phi_ip_km_lm = Root_Mesh_Hydro[index_root(i,k,l)].gravity
								[Root_Mesh_Hydro[index_root(i,k,l)].
								 index_nested_gravity_mesh(i_nst+1,k_nst  ,l_nst  )];
							phi_im_kp_lp = Root_Mesh_Hydro[index_root(i,k,l)].gravity
								[Root_Mesh_Hydro[index_root(i,k,l)].
								 index_nested_gravity_mesh(i_nst  ,k_nst+1,l_nst+1)];
							phi_im_kp_lm = Root_Mesh_Hydro[index_root(i,k,l)].gravity
								[Root_Mesh_Hydro[index_root(i,k,l)].
								 index_nested_gravity_mesh(i_nst  ,k_nst+1,l_nst  )];
							phi_im_km_lp = Root_Mesh_Hydro[index_root(i,k,l)].gravity
								[Root_Mesh_Hydro[index_root(i,k,l)].
								 index_nested_gravity_mesh(i_nst  ,k_nst  ,l_nst+1)];
							phi_im_km_lm = Root_Mesh_Hydro[index_root(i,k,l)].gravity
								[Root_Mesh_Hydro[index_root(i,k,l)].
								 index_nested_gravity_mesh(i_nst  ,k_nst  ,l_nst  )];

							// ���������� ���������
							gradfix = (
								0.25*(	phi_ip_kp_lp + phi_ip_kp_lm + 
										phi_ip_km_lp + phi_ip_km_lm) - 
								0.25*(	phi_im_kp_lp + phi_im_kp_lm + 
										phi_im_km_lp + phi_im_km_lm) )/h_nst;
							gradfiy = (
								0.25*(	phi_ip_kp_lp + phi_ip_kp_lm +
										phi_im_kp_lp + phi_im_kp_lm) -
								0.25*(	phi_ip_km_lp + phi_ip_km_lm +
										phi_im_km_lp + phi_im_km_lm) )/h_nst;
							gradfiz = (
								0.25*(	phi_ip_kp_lp + phi_ip_km_lp + 
										phi_im_kp_lp + phi_im_km_lp) -
								0.25*(	phi_ip_kp_lm + phi_ip_km_lm + 
										phi_im_kp_lm + phi_im_km_lm) )/h_nst;

							// �������������� ��������
							cell.rhovx = cell.rhovx - tau * cell.rho * gradfix;
							cell.rhovy = cell.rhovy - tau * cell.rho * gradfiy;
							cell.rhovz = cell.rhovz - tau * cell.rho * gradfiz;

							// ������ ��������� �������� � ������
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhovx = 
									cell.rhovx;
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhovy = 
									cell.rhovy;
							Root_Mesh_Hydro[index_root(i,k,l)].
								mesh[Root_Mesh_Hydro[index_root(i,k,l)].
								index_nested_mesh(i_nst,k_nst,l_nst)].rhovz = 
									cell.rhovz;
						}
			}
	
}
/*******************************************************************/
#endif