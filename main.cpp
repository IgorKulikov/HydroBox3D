/* Проект HydroBox3D */
#include "system.hpp"
#include "global.hpp"
#include "setup.hpp"
#include "output.hpp"
#include "mesh.hpp"
#include "cfl.hpp"
#include "boundary.hpp"
#include "godunov.hpp"
#include "poisson.hpp"
#include "physics.hpp"

int main()
{
	double timer = 0.0, tend;

	Init();						// создание массивов и окружения
	Setup();					// формулировка задачи
	Primitive_Variable();		// вычисление физических переменных
	Poisson_Solver();			// решение уравнения Пуассона
	Save(timer);				// сохранение результата
	LogWrite(timer);			// сохранение начальных значений энергий
	
	// моделирование гидродинамического процесса
	for(int i=0 ; i<iTimeCount ; i++)
	{
		tend = print_tau*(i+1);		// вычисление времени подэтапа
		Mesh_Reconstruct();			// перестройка сетки
		while(timer < tend)
		{
			Primitive_Variable();		// вычисление физических переменных
			Poisson_Solver();			// решение уравнения Пуассона
			Tau_Computing(timer,tend);	// вычисление шага tau
			Boundary_Condition();		// постановка граничных условий
			Riemann_Computing();		// вычисление потоков Римана
			Godunov_Solver();			// реализация схемы Годунова
			Gravity_Force();			// вычисление гравитационной силы
			timer += tau;
			LogWrite(timer);			// сохранение значений энергий
			printf("*",timer);
		}
		printf("\n");
		Save(timer);				// сохранение результата
	}

	Destroy();		// удаление массивов и окружения
	return 0;
}