// Individual.h: interface for the CIndividual class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_INDIVIDUAL_H__A05E59DB_40D8_4282_9D3C_A0201898ACA8__INCLUDED_)
#define AFX_INDIVIDUAL_H__A05E59DB_40D8_4282_9D3C_A0201898ACA8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "../StdAfx.h"

/* the following parameters are problem dependent */
#define index_of_normal			1					// index of the normal function, 0 means no normal function


#define N_of_x					100					// number of the decision variables
#define Max_of_NFEs				(N_of_x*5000)		// the maximal number of function evaluations

// for constrained problems
#define N_of_constr				0					// number of the constrained functions

#define PRECISION				1e-8				// value to reach (VTR)
/* ---------------------------------------------- */

#define population_size			(4*N_of_x)			// size of the population
/* -------------------------------------------------------------- */

class CIndividual  
{
public:
	CIndividual();
	virtual ~CIndividual();

public:
	CIndividual(const CIndividual &);							// 拷贝构造函数
	CIndividual &operator=(const CIndividual &);				// 重载"="运算符

public:
	double		*xreal;				// Real Coded Decision Variable     
	tFitness	obj;				// value of the objective function
	tFitness	fitness;			// fitness of the objective function
	double		constr_violation;	// parameter for constraint violation
    double		*constr;			// defining the constraint values
	int			feasible;			// check the individual is feasible or not (1: feasible; 0: infeasible)
	int			no_of_violation;	// number of the violated constraint functions

	double		F;					// only for DE
	double		CR;					// only for DE
};

#endif // !defined(AFX_INDIVIDUAL_H__A05E59DB_40D8_4282_9D3C_A0201898ACA8__INCLUDED_)
