// ProblemDef.h: interface for the CProblemDef class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PROBLEMDEF_H__C8A8E61F_A57B_4F2A_A52A_0956308C172C__INCLUDED_)
#define AFX_PROBLEMDEF_H__C8A8E61F_A57B_4F2A_A52A_0956308C172C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Individual.h"
#include "Rand.h"


// only for unconstrained optimization problems
class CProblemDef  
{
public:
	CProblemDef();
	virtual ~CProblemDef();

	CRand m_rnd;

	// for unconstrained optimization problems
public:
	double	center[200];
	int		flag_once;
	double	TempValue(double x,int a,int k,int m);

	void	evaluate_normal_fitness(double *xreal,tFitness &obj, 
		double *constr, int func_flag, long int &evaluations);

public:
	void test_01(double *xreal,tFitness &obj, double *constr);
	void test_02(double *xreal,tFitness &obj, double *constr);
	void test_03(double *xreal,tFitness &obj, double *constr);
	void test_04(double *xreal,tFitness &obj, double *constr);
	void test_05(double *xreal,tFitness &obj, double *constr);
	void test_06(double *xreal,tFitness &obj, double *constr);
	void test_07(double *xreal,tFitness &obj, double *constr);
	void test_08(double *xreal,tFitness &obj, double *constr);
	void test_09(double *xreal,tFitness &obj, double *constr);
	void test_10(double *xreal,tFitness &obj, double *constr);
	void test_11(double *xreal,tFitness &obj, double *constr);
	void test_12(double *xreal,tFitness &obj, double *constr);
	void test_13(double *xreal,tFitness &obj, double *constr);
	
};

#endif // !defined(AFX_PROBLEMDEF_H__C8A8E61F_A57B_4F2A_A52A_0956308C172C__INCLUDED_)
