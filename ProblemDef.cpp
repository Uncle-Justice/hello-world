// ProblemDef.cpp: implementation of the CProblemDef class.
//
//////////////////////////////////////////////////////////////////////

// #include "stdafx.h"
#include "ProblemDef.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CProblemDef::CProblemDef()
{
}

CProblemDef::~CProblemDef()
{
}

void CProblemDef::evaluate_normal_fitness(double *xreal, tFitness &obj, double *constr, int func_flag, long int &evaluations)
{
	switch(func_flag) {
	case 1:
		test_01(xreal,obj,constr);
		break;
	case 2:
		test_02(xreal,obj,constr);
		break;
	case 3:
		test_03(xreal,obj,constr);
		break;
	case 4:
		test_04(xreal,obj,constr);
		break;
	case 5:
		test_05(xreal,obj,constr);
		break;
	case 6:
		test_06(xreal,obj,constr);
		break;
	case 7:
		test_07(xreal,obj,constr);
		break;
	case 8:
		test_08(xreal,obj,constr);
		break;
	case 9:
		test_09(xreal,obj,constr);
		break;
	case 10:
		test_10(xreal,obj,constr);
		break;
	case 11:
		test_11(xreal,obj,constr);
		break;
	case 12:
		test_12(xreal,obj,constr);
		break;
	case 13:
		test_13(xreal,obj,constr);
		break;	
	default:
		printf("The function you selected does not exist.\n");
		exit(0);
		break;
	}
	evaluations++;
}

// only for unconstrained optimization problems
// Ellipsoidal function
void CProblemDef::test_01(double *x, tFitness &obj, double *constr)
{
	tFitness fit = 0.0;
	for(int i=0;i<N_of_x;i++)
		fit += ((x[i])*(x[i]));
	obj = fit;
}

// Schwefel's function 2.22
void CProblemDef::test_02(double *x, tFitness &obj, double *constr)
{
	tFitness value = 0.0;
	tFitness temp1 = 0.0;
	tFitness temp2 = 1.0;
	for (int i=0;i<N_of_x;i++)
	{
		 temp1 += fabs(x[i]);
		 temp2 *= fabs(x[i]);
	}
	value = temp1+temp2;
	obj = value;
}

// Schwefel's function 1.2
void CProblemDef::test_03(double *x, tFitness &obj, double *constr)
{
	tFitness fit = 0.0;
	tFitness sumSCH;
	int i,j;
	
	for (i=0;i<N_of_x;i++)
	{
		sumSCH = 0.0;
		for (j=0;j<i+1;j++)
		{
			sumSCH += x[j];
		}
		fit += sumSCH*sumSCH;
	}
	obj = fit;
}

// Schwefel's function 2.21
void CProblemDef::test_04(double *x, tFitness &obj, double *constr)
{
	tFitness temp1 = fabs(x[0]);
	for (int i=1;i<N_of_x;i++)
	{
		if (fabs(x[i]) > temp1)
		{
			temp1 = fabs(x[i]);
		}
	}
	obj = temp1;
}

// Rosenbrock's function
void CProblemDef::test_05(double *x, tFitness &obj, double *constr)
{
	tFitness fit;
	tFitness t0, tt, t1, d=0;
	t0=x[0];
	for (int i=1; i<N_of_x; i++) 
	{
		t1 = x[i];
		tt = (1.0-t0);
		d += tt*tt;
		tt = t1-t0*t0;
		d += 100*tt*tt;				
		t0 = t1;
	}
	fit = d;

	obj = fit;
}

// step function
void CProblemDef::test_06(double *x, tFitness &obj, double *constr)
{
	tFitness temp;
	tFitness value = 0.0;
	for (int i=0;i<N_of_x;i++)
	{
		temp = floor(x[i]+0.5);
		value += temp*temp; 
	}
	obj = value;
}

// Quartic function
void CProblemDef::test_07(double *xreal, tFitness &obj, double *constr)
{
	tFitness value = 0.0;
	for (int i=0;i<N_of_x;i++)
	{
		value += (i+1)*pow(xreal[i],4.0); 
	}
	CRand rnd;
	obj = value+rnd.randomperc();
}

// Generalized Schwefel's function
void CProblemDef::test_08(double *x, tFitness &obj, double *constr)
{
	tFitness value = 0.0;
	for (int i=0;i<N_of_x;i++)
	{
		value += x[i]*sin(sqrt(fabs(x[i]))); 
	}
	obj = -value + 418.9828872724337998*(double)N_of_x; // modified here
}

// Generalized Rastrigin's function
void CProblemDef::test_09(double *x, tFitness &obj, double *constr)
{
	tFitness value = 0.0;
	for (int i=0;i<N_of_x;i++)
	{
		value += x[i]*x[i]-10.0*cos(2.0*PI*x[i])+10.0; 
	}
	obj = value;
}

// Ackley' function
void CProblemDef::test_10(double *x, tFitness &obj, double *constr)
{
	tFitness value = 0.0;
	tFitness temp1 = 0.0;
	tFitness temp2 = 0.0;
	for (int i=0;i<N_of_x;i++)
	{
		 temp1 += x[i]*x[i];
		 temp2 += cos(2.0*PI*x[i]);
	}
	obj = -20.0*exp(-0.2*sqrt(temp1/N_of_x))-exp(temp2/N_of_x)+20.0+exp(1.0);
}

// Generalized Griewank function
void CProblemDef::test_11(double *x, tFitness &obj, double *constr)
{
	tFitness temp1 = 0.0;
	tFitness temp2 = 1.0;
	for (int i=0;i<N_of_x;i++)
	{
		 temp1 += x[i]*x[i];
		 temp2 *= cos(x[i]/sqrt(i+1));
	}
	obj = temp1/4000.0-temp2+1;
}

double CProblemDef::TempValue(double x,int a,int k,int m)
{
	double temp = 0.0;
	if( x > a)
	{
		temp = k*pow(x-a,m);
	}
	else if( x <= a && x >= -a)
	{
		temp = 0.0;
	}
	else
	{
		temp = k*pow(-x-a,m);
	}
	return temp;
}

// Generalized Penalized function
void CProblemDef::test_12(double *x, tFitness &obj, double *constr)
{
	double *y;	//��ʱ�洢�����еı���

	y=new double[N_of_x];

	for (int i=0;i<N_of_x;i++)
	{
		y[i]=0.0;
	}

	for (int i=0;i<N_of_x;i++)
	{
		y[i]=1+(x[i]+1)/4.0;
	}

	tFitness temp1 = 0.0;
	tFitness temp2 = 0.0;
	for (int i=0;i<N_of_x-1;i++)
	{
		temp1 += pow(y[i]-1,2.0)*(1.0+10.0*pow(sin(PI*y[i+1]),2.0)); 
	}
	for (int i=0;i<N_of_x;i++)
	{
		temp2 += TempValue(x[i],10,100,4);
	}
	obj = (10.0*pow(sin(PI*y[0]),2.0)+temp1+pow(y[N_of_x-1]-1,2))*PI/N_of_x+temp2;
	delete []y;
}

// Generalized Penalized function
void CProblemDef::test_13(double *x, tFitness &obj, double *constr)
{
	tFitness temp1 = 0.0;
	tFitness temp2 = 0.0;
	for (int i=0;i<N_of_x-1;i++)
	{
		temp1 += pow(x[i]-1,2.0)*(1.0+10.0*pow(sin(3*PI*x[i+1]),2.0)); 
	}
	for (int i=0;i<N_of_x;i++)
	{
		temp2 += TempValue(x[i],5,100,4);
	}
	obj = (pow(sin(3.0*PI*x[0]),2.0)+temp1+pow(x[N_of_x-1]-1,2.0)
		*(1.0+pow(sin(2.0*PI*x[N_of_x-1]),2.0)))/10.0+temp2;
}


