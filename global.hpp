#ifndef _global_hpp_
#define _global_hpp_
/*******************************************************************/
#include "system.hpp"

// ����� ��
#define PI 3.141592653589793

// ��������� �������� �����
#define N 32			// ������ �������� �����
#define PSN_POW 5		// ������� ������ �������� ����� [N = 2^PSN_POW]
#define MAX_NESTED_N 8	// ������������ ������ ��������� �����

// ������������ ����� �������� ������� ��������� �������� ������� SOR
#define MAX_SOR_ITER 10

// �������� ���������� ������ SOR
#define RELAXATION_SOR 0.7

const double box_size = 3.2;			// ������ �������
double root_h = (box_size/(double)N);	// ��� �������� �����

const double print_tau = 0.1;	// ��� ����������� �����
const int iTimeCount = 30;		// ����� ����� ��������

int root_nx = N;		// ��������� ���������� �������� �����
int root_ny = N;
int root_nz = N;

// ��� �� ������� � ������� �������
#define CFL 0.2
double tau;

// ��������� �������
struct root_cell *Root_Mesh_Hydro;

// ��������������� ������� ��� ������� ��������� ��������
double *xfft, *yfft;

// ���� ��� �������� �� ���������� �������
FILE *fout;

/*******************************************************************/
#endif