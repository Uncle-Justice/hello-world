// Individual.cpp: implementation of the CIndividual class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
#include "Individual.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CIndividual::CIndividual()
{
	xreal = new double[N_of_x];
	constr = new double[N_of_constr];
	for (int i=0;i<N_of_x;i++)
	{
		xreal[i] = -1.0;
	}
	for (int i=0;i<N_of_constr;i++)
	{
		constr[i] = -1.0;
	}
	constr_violation = -1.0;
	obj              = -1.0;
	fitness          = -1.0;
	feasible         = 0;
	no_of_violation  = 0;

	F  = 0.5;
	CR = 0.9;
}

CIndividual::~CIndividual()
{
	delete []xreal;
	if((N_of_constr != 0)) delete[] constr;
}

CIndividual::CIndividual(const CIndividual &t)
{
	// �ж��Ƿ���ͬ
	if(&t == this)
	{
		return;
	}
	
	// ��ֵ��this
	xreal = new double[N_of_x];
	constr = new double[N_of_constr];
	for (int i=0;i<N_of_x;i++)
	{
		xreal[i] = t.xreal[i];
	}
	obj = t.obj;
	fitness = t.fitness;
	for (int i=0;i<N_of_constr;i++)
	{
		constr[i] = t.constr[i];
	}
	constr_violation = t.constr_violation;
	feasible         = t.feasible;
	no_of_violation  = t.no_of_violation;

	F  = t.F;
	CR = t.CR;
}

CIndividual &CIndividual::operator =(const CIndividual &t)
{
	// �ж��Ƿ���ͬ
	if(&t == this)
	{
		return *this;
	}

    // !!!!�ж��Ƿ��Ѿ�����ռ�,����Ҫ,��Ȼ�ڴ治�������ͷ�!!!!
	if ((N_of_x != 0) && (xreal)) delete[] xreal;
	if((N_of_constr != 0) && (constr)) delete[] constr;

	xreal = new double[N_of_x];
	constr = new double[N_of_constr];
	for (int i=0;i<N_of_x;i++)
	{
		xreal[i] = t.xreal[i];
	}
	obj = t.obj;
	fitness = t.fitness;
	for (int i=0;i<N_of_constr;i++)
	{
		constr[i] = t.constr[i];
	}
	constr_violation = t.constr_violation;
	feasible         = t.feasible;
	no_of_violation  = t.no_of_violation;

	F  = t.F;
	CR = t.CR;

	return *this;
}
