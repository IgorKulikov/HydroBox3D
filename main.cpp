/* ������ HydroBox3D */
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

	Init();						// �������� �������� � ���������
	Setup();					// ������������ ������
	Primitive_Variable();		// ���������� ���������� ����������
	Poisson_Solver();			// ������� ��������� ��������
	Save(timer);				// ���������� ����������
	LogWrite(timer);			// ���������� ��������� �������� �������
	
	// ������������� ������������������ ��������
	for(int i=0 ; i<iTimeCount ; i++)
	{
		tend = print_tau*(i+1);		// ���������� ������� ��������
		Mesh_Reconstruct();			// ����������� �����
		while(timer < tend)
		{
			Primitive_Variable();		// ���������� ���������� ����������
			Poisson_Solver();			// ������� ��������� ��������
			Tau_Computing(timer,tend);	// ���������� ���� tau
			Boundary_Condition();		// ���������� ��������� �������
			Riemann_Computing();		// ���������� ������� ������
			Godunov_Solver();			// ���������� ����� ��������
			Gravity_Force();			// ���������� �������������� ����
			timer += tau;
			LogWrite(timer);			// ���������� �������� �������
			printf("*",timer);
		}
		printf("\n");
		Save(timer);				// ���������� ����������
	}

	Destroy();		// �������� �������� � ���������
	return 0;
}