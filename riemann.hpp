#ifndef _riemann_hpp_
#define _riemann_hpp_
/*******************************************************************/
#include <math.h>
#include "system.hpp"
#include "physics.hpp"

// римановский решатель
flux riemann(double rho_minus,  double rho_plus, 
			 double velx_minus, double velx_plus, 
			 double vely_minus, double vely_plus, 
			 double velz_minus, double velz_plus, 
			 double p_minus,    double p_plus)
{
	double density, pressure, velx, vely, velz, velocity, csound, maxlambda, 
		   ent_minus, ent_plus, entropy;
	double frho_minus, frho_plus, frhovx_minus, frhovx_plus, 
		   frhovy_minus, frhovy_plus, frhovz_minus, frhovz_plus,
		   frhoent_minus, frhoent_plus;
	flux riemann_solution;

	// перевычисление энтропии
	ent_minus = get_entropy(rho_minus,p_minus);
	ent_plus  = get_entropy(rho_plus ,p_plus );

	// вычисление значений на интерфейсе
	density  = (sqrt(rho_minus)*rho_minus+sqrt(rho_plus)*rho_plus) / 
				(sqrt(rho_minus)+sqrt(rho_plus));
	pressure = (sqrt(rho_minus)*p_minus+sqrt(rho_plus)*p_plus) / 
				(sqrt(rho_minus)+sqrt(rho_plus));
	velx     = (sqrt(rho_minus)*velx_minus+sqrt(rho_plus)*velx_plus) / 
				(sqrt(rho_minus)+sqrt(rho_plus));
	vely     = (sqrt(rho_minus)*vely_minus+sqrt(rho_plus)*vely_plus) / 
				(sqrt(rho_minus)+sqrt(rho_plus));
	velz     = (sqrt(rho_minus)*velz_minus+sqrt(rho_plus)*velz_plus) / 
				(sqrt(rho_minus)+sqrt(rho_plus));
	entropy  = (sqrt(rho_minus)*ent_minus+sqrt(rho_plus)*ent_plus) / 
				(sqrt(rho_minus)+sqrt(rho_plus));
	
	velocity  = sqrt(velx*velx+vely*vely+velz*velz);
	csound    = sqrt(ADIABATIC_INDEX_GAMMA*pressure/density);
	maxlambda = csound + velocity;
		
	// формирование потоков
	frho_minus   = rho_minus * velx_minus;
	frho_plus    = rho_plus  * velx_plus;
	
	frhovx_minus = rho_minus * velx_minus * velx_minus + p_minus;
	frhovx_plus  = rho_plus  * velx_plus  * velx_plus  + p_plus;
	frhovy_minus = rho_minus * vely_minus * velx_minus;
	frhovy_plus  = rho_plus  * vely_plus  * velx_plus;
	frhovz_minus = rho_minus * velz_minus * velx_minus;
	frhovz_plus  = rho_plus  * velz_plus  * velx_plus;

	frhoent_minus = rho_minus * ent_minus * velx_minus;
	frhoent_plus  = rho_plus  * ent_plus  * velx_plus;

	// HLLK решатель
	riemann_solution.flux_rho   = 0.5*(frho_minus + frho_plus) + 
		0.5*maxlambda*(rho_minus - rho_plus);
	riemann_solution.flux_rhovx = 0.5*(frhovx_minus + frhovx_plus) + 
		0.5*maxlambda*(rho_minus*velx_minus - rho_plus*velx_plus);
	riemann_solution.flux_rhovy = 0.5*(frhovy_minus + frhovy_plus) + 
		0.5*maxlambda*(rho_minus*vely_minus - rho_plus*vely_plus);
	riemann_solution.flux_rhovz = 0.5*(frhovz_minus + frhovz_plus) + 
		0.5*maxlambda*(rho_minus*velz_minus - rho_plus*velz_plus);
	riemann_solution.flux_rhoent = 0.5*(frhoent_minus + frhoent_plus) + 
		0.5*maxlambda*(rho_minus*ent_minus - rho_plus*ent_plus);

	return riemann_solution;
}
/*******************************************************************/
#endif