#ifndef _global_hpp_
#define _global_hpp_
/*******************************************************************/
#include "system.hpp"

// Число Пи
#define PI 3.141592653589793

// Параметры корневой сетки
#define N 32			// размер корневой сетки
#define PSN_POW 5		// степень двойки корневой сетки [N = 2^PSN_POW]
#define MAX_NESTED_N 8	// максимальный размер вложенной сетки

// Максимальное число итераций решения уравнения Пуассона методом SOR
#define MAX_SOR_ITER 10

// Параметр релаксации метода SOR
#define RELAXATION_SOR 0.7

const double box_size = 3.2;			// размер области
double root_h = (box_size/(double)N);	// шаг корневой сетки

const double print_tau = 0.1;	// шаг перестройки сетки
const int iTimeCount = 30;		// число шагов процесса

int root_nx = N;		// настройка параметров корневой сетки
int root_ny = N;
int root_nz = N;

// Шаг по времени и условие Куранта
#define CFL 0.2
double tau;

// Расчетная область
struct root_cell *Root_Mesh_Hydro;

// Вспомогательные массивы для решения уравнения Пуассона
double *xfft, *yfft;

// Файл для контроля за поведением энергий
FILE *fout;

/*******************************************************************/
#endif