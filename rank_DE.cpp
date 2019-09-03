// DE_probability.cpp : Defines the entry point for the console application.
//
#pragma warning(disable:4996)
#include "StdAfx.h"
#include ".\DE_method.h"
using namespace std;
void welcome();

int main(int argc, char* argv[])
{
	welcome();

	double seed = 0.0;

	int method;
	int RUN_NUMBER;
	int strategy_index = 1;		// strategy used in DE (1: rand/1; )
	int output_counter = 0;
	
	/**/printf("method = ");
	method = 1;
	cin>>method;

	printf("\nPlease enter the max run number.\n");
	printf("RUN_NUMBER = ");
	RUN_NUMBER = 5;
	cin>>RUN_NUMBER;
	
	/*if (argc != 4 && argc != 5)
	{
		printf("Usage: DE_*.exe method run_number out_counter\n");
		printf("OR--->\n");
		printf("Usage: DE_*.exe method run_number out_counter strategy_index\n");
		exit(0);
	}
	if (argc == 4)
	{
		method         = atoi(argv[1]);
		RUN_NUMBER     = atoi(argv[2]);
		output_counter = atoi(argv[3]);
	}
	else if (argc == 5)
	{
		method         = atoi(argv[1]);
		RUN_NUMBER     = atoi(argv[2]);
		output_counter = atoi(argv[3]);
		strategy_index = atoi(argv[4]);
 	}*/

	srand((unsigned)time(NULL));


	CDE_method *DE_method;
	for (int i=0;i<RUN_NUMBER;i++)
	{
		int t = output_counter;				// ��������̫����������ô����з���
		int k = t*RUN_NUMBER;
		printf("--------- Run no. is %d ---------\n", k+i+1);
		//seed = ((double)(k+i+1))/((double)RUN_NUMBER);
		seed = ((double)(k+i+1))/50.0;		// ............
		if (RUN_NUMBER == 1)
		{
			//seed = 0.86564;
			seed = 0.2;
		}
		//seed = ((double)(46+1))/50.0;			// test seed
		//seed = (double)(rand()%1000)/1000.0;	// test seed
		if (seed == 0.0)
		{
			seed = 0.001;
		}
		if (seed == 1.0)
		{
			seed = 0.99;
		}

		DE_method = new CDE_method;
		DE_method->Run_Optimizer(k+i, seed, method, strategy_index);
		delete DE_method;

		if ((i+1)%100 == 0 && (i+1) != RUN_NUMBER)	// pause
		{
			printf("\nPlease press ENTER key to continue.\n");
			getchar();
		}
	}
	
	return 0;
}

void welcome()
{
	printf("\n");
	printf("***************************************************************\n");
	printf("*             Welcome to use the DE optimizers.               *\n");
	printf("*        You can select and use the following optimizer.      *\n");
	printf("*                method = 1   for DE                          *\n");
	printf("*                method = 2   for rank_DE                     *\n");
	printf("***************************************************************\n");
	printf("\n");
}
