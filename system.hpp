#ifndef _system_hpp_
#define _system_hpp_
/*******************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.hpp"

/* ������� ������� �������� ������ */
int index_root(int i, int k, int l)
{
	return i*root_nz*root_ny+k*root_nz+l;
}

/* ��������� ������ ����� ������� */
struct flux
{
	double flux_rho, flux_rhovx, flux_rhovy, flux_rhovz, flux_rhoent;
};

/* ��������� ������ ��������� ����� */
struct nested_cell	
{
	double rho,						// ���������� ����������
		   velx, vely, velz, 
		   entropy,
		   press, 
		   rhoent, rhovx, rhovy, rhovz;		// �������������� ����������
	
	// ������ ����� �������
	flux   x_plus, x_minus, y_plus, y_minus, z_plus, z_minus;
};

/* ��������� ������ �������� ������ */
struct root_cell		
{
	double rho,			// ������� ��������� - �������� ��������� ����� � ����������
		   phi,			// �������������� ���������
		   phi_real,	// ������� ��� �������������� �����
		   phi_image;

	int nested_N;		// ������ ��������� �����
	double nested_h;	// ��� ��������� ����� [nested_h = root_h/nested_N]
	double *gravity;	// ������� ��������� ����� ��� ���������� (nested_N+1)^3
	double *grav_new;	// ������� ��������� ����� ��� ������������ ��������
	nested_cell *mesh;	// ��������� ����� ��� ��������� ������������� (nested_N)^3
	
	int index_nested_gravity_mesh(int i, int k, int l)	// ������ ������ ����������
	{
		return i*(nested_N+1)*(nested_N+1)+k*(nested_N+1)+l;
	}

	int index_nested_mesh(int i, int k, int l)	// ������ ������ ��������� �����
	{
		return i*nested_N*nested_N+k*nested_N+l;
	}

	void make_nested_mesh(int new_N)	// �������� ��������� �����
	{
		if(new_N > MAX_NESTED_N) new_N = MAX_NESTED_N;
		nested_N = new_N;
		nested_h = root_h/nested_N;
		mesh = new nested_cell[nested_N*nested_N*nested_N];
		gravity = new double[(nested_N+1)*(nested_N+1)*(nested_N+1)];
		grav_new = new double[(nested_N+1)*(nested_N+1)*(nested_N+1)];
	}

	void reconstruct_nested_mesh(int new_N)	// ������������� ��������� �����
	{
		return;		/*** �������-�������� ***/
	}

	void destroy_nested_mesh()	// �������� ��������� �����
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