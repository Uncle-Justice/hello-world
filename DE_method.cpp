// DE_method.cpp: implementation of the CDE_method class.
//
//////////////////////////////////////////////////////////////////////
#pragma warning(disable:4996)
//#include "stdafx.h"
#include ".\DE_method.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CDE_method::CDE_method()
{
	max_iterations = Max_of_NFEs/population_size;
	max_NFFEs      = Max_of_NFEs;
	pop_size       = population_size;

	parent_pop   = new CIndividual[2*pop_size];
	child_pop    = new CIndividual[2*pop_size];
	child_pop1   = new CIndividual[2*pop_size];
	upper_bounds = new double[N_of_x];
	lower_bounds = new double[N_of_x];
	sort_index   = new int[2*pop_size];

	/************************************************************************/
	/* default parameter settings for DE method                             */
	/************************************************************************/
	m_F              = 0.5;
	m_CR             = 0.9;
	m_strategy       = 1;		// DE/rand/1/bin
	param_adaptation = 1;		// for parameter adaptation 

	/************************************************************************/
	/* for probability-based DE methods                                     */
	/************************************************************************/
	weights_flag     = 2;			// 1: ranking-based linear model; 2: ranking-based quadratic model; 3: ranking-based cosine model
}

CDE_method::~CDE_method()
{
	delete []parent_pop;
	delete []child_pop;
	delete []child_pop1;
	delete []upper_bounds;
	delete []lower_bounds;
	delete []sort_index;
}

void CDE_method::init_variables()
{
	if (index_of_normal != 0)
	{
		init_normal_variables();
	}

	max_interval = 0;
	double temp;
	for (int i=0;i<N_of_x;i++)
	{
		temp = upper_bounds[i]-lower_bounds[i];
		max_interval += temp*temp;
	}
	max_interval = sqrt(max_interval);
}

void CDE_method::init_normal_variables()
{
	int i;
	int n     = N_of_x;
	func_flag = index_of_normal;

	switch(func_flag) {
	case 1:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -100.0;
			upper_bounds[i] = 100.0;
		}
		known_optimal = 0.0;
		MINIMIZE      = 1;
		break;
	case 2:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -10.0;
			upper_bounds[i] = 10.0;
		}
		known_optimal = 0.0;
		MINIMIZE      = 1;
		break;
	case 3:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -100.0;
			upper_bounds[i] = 100.0;
		}
		known_optimal = 0.0;
		MINIMIZE      = 1;
		break;
	case 4:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -100.0;
			upper_bounds[i] = 100.0;
		}
		known_optimal = 0.0;
		MINIMIZE      = 1;
		break;
	case 5:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -30.0;
			upper_bounds[i] = 30.0;
		}
		known_optimal = 0.0;
		MINIMIZE      = 1;
		break;
	case 6:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -100.0;
			upper_bounds[i] = 100.0;
		}
		known_optimal = 0.0;
		MINIMIZE      = 1;
		break;
	case 7:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -1.28;
			upper_bounds[i] = 1.28;
		}
		known_optimal = 0.0;
		MINIMIZE      = 1;
		break;
	case 8:// -418.982887272434*D
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -500.0;
			upper_bounds[i] = 500.0;
		}
		known_optimal = 0.0;	// 7.27595761418343e-012
		MINIMIZE      = 1;
		break;
	case 9:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -5.12;
			upper_bounds[i] = 5.12;
		}
		known_optimal = 0.0;
		MINIMIZE      = 1;
		break;
	case 10:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -32.0;
			upper_bounds[i] = 32.0;
		}
		known_optimal = 5.88721779659630e-016;
		MINIMIZE      = 1;
		break;
	case 11:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -600.0;
			upper_bounds[i] = 600.0;
		}
		known_optimal = 0.0;
		MINIMIZE      = 1;
		break;
	case 12:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -50.0;
			upper_bounds[i] = 50.0;
		}
		known_optimal = 1.57044103551778e-032;
		MINIMIZE      = 1;
		break;
	case 13:
		for (i=0;i<N_of_x;i++)
		{
			lower_bounds[i] = -50.0;
			upper_bounds[i] = 50.0;
		}
		known_optimal = 1.34969464963993e-032;
		MINIMIZE      = 1;
		break;
	default:
		printf("The function you selected does not exist.\n");
		exit(0);
	}
}

void CDE_method::init_pop()
{
	printf("Random based parent_pop initialization method is used.\n");
	init_pop_random();

	// initialize the F and CR of DE
	for (int i=0;i<pop_size;i++)
	{
		parent_pop[i].F  = m_rnd.rndreal(0.1, 1.0);
		parent_pop[i].CR = m_rnd.rndreal(0.0, 1.0);
	}
}

void CDE_method::init_pop_random()
{
	int i, j;
	for (i=0;i<pop_size;i++)
	{
		// for bounded problem
		for (j=0;j<N_of_x;j++)
		{
			parent_pop[i].xreal[j] = m_rnd.rndreal(lower_bounds[j],upper_bounds[j]);
		}
	}
}

void CDE_method::evaluate_ind(CIndividual &indv)
{
	if (index_of_normal != 0)
	{
		evaluate_normal_ind(indv);
	}

	// convert objective function value into fitness
	indv.fitness = MINIMIZE*indv.obj;

	// sum violation of the constrained functions
	if (N_of_constr == 0)
	{
		indv.constr_violation = 0.0;
		indv.feasible = 1;
	}
}

void CDE_method::evaluate_normal_ind(CIndividual &indv)
{
	func_flag = index_of_normal;

	m_func.evaluate_normal_fitness(indv.xreal,indv.obj,indv.constr, func_flag, evaluations);
}

void CDE_method::evaluate_pop(CIndividual *pop, int size)
{
	for (int i=0;i<size;i++)
	{
		evaluate_ind(pop[i]);
	}
}

int CDE_method::compare_ind (CIndividual *indv1, CIndividual *indv2)
{
	if((indv1->feasible==TRUE && indv2->feasible==TRUE))
	{
		if(indv1->fitness < indv2->fitness)
		{
			return 1;
		}
		else if (indv1->fitness > indv2->fitness)
		{
			return -1;
		}
		else
		{
			return 0;
		}
	}
	else if(indv1->feasible==TRUE && indv2->feasible==FALSE)
	{
		return 1;
	}
	else if (indv1->feasible==FALSE && indv2->feasible==TRUE)
	{
		return -1;
	}
	else
	{
		if(indv1->constr_violation < indv2->constr_violation)
		{
			return 1;
		}
		else if (indv1->constr_violation > indv2->constr_violation)
		{
			return -1;
		}
		else
		{
			//return 0;
			if (indv1->fitness < indv2->fitness)
			{
				return 1;
			}
			else if (indv1->fitness > indv2->fitness)
			{
				return -1;
			}
			else
			{
				return 0;
			}
		}
	}
}

void CDE_method::find_best_index(CIndividual *pop, int size)
{
	int flag;

	best_index = 0;
	for (int i=1;i<size;i++)
	{
		flag = compare_ind(&pop[i], &pop[best_index]);
		if (flag == 1)
		{
			best_index = i;
		}
	}

	worst_index = 0;
	for (int i =1;i<size;i++)
	{
		flag = compare_ind(&pop[i], &pop[worst_index]);
		if (flag == -1)
		{
			worst_index = i;
		}
	}
}

void CDE_method::random_index(int *array_index, int all_size, int size)
{
	int i,j,krand;
	int *a;
	a = new int[all_size];

	for(i = 0;i<all_size;i++)
	{
		a[i] = i;
	}
	for(i=0;i<size;i++)
	{
		j     = m_rnd.rndint(i,(all_size-1));
		krand = a[i];
		a[i]  = a[j];
		a[j]  = krand;
	}	
	for(i=0;i<size;i++)
	{
		array_index[i] = a[i];
	}

	//a = NULL;
	delete []a;
}

void CDE_method::shell_sort_pop(CIndividual *pop, int *list, int size)
{
	int done;
	int step, bound, i, j;
	int temp;

	for (i=0;i<size;i++)
	{
		list[i] = i;
	}

	step = size;  // array length
	while (step > 1) 
	{
		step /= 2;	//halve the step size
		do 
		{
			done   = 1;
			bound  = size - step;
			for (j = 0; j < bound; j++) 
			{
				i = j + step + 1;
				if (compare_ind(&pop[list[j]], &pop[list[i-1]]) == -1) 	
				{
					temp      = list[i-1];
					list[i-1] = list[j];
					list[j]   = temp;
					done      = 0; // if a swap has been made we are not finished yet
				}  // if
			}  // for
		} while (done == 0);   // while
	} //while (step > 1)
}

void CDE_method::shell_sort_array(double *array, int *list, int size)
{
	int done;
	int step, bound, i, j;
	int temp;

	for (i=0;i<size;i++)
	{
		list[i] = i;
	}

	step = size;  // array length
	while (step > 1) 
	{
		step /= 2;	//halve the step size
		do 
		{
			done   = 1;
			bound  = size - step;
			for (j = 0; j < bound; j++) 
			{
				i = j + step + 1;
				if (array[list[j]] < array[list[i-1]]) 	
				{
					temp      = list[i-1];
					list[i-1] = list[j];
					list[j]   = temp;
					done      = 0; // if a swap has been made we are not finished yet
				}  // if
			}  // for
		} while (done == 0);   // while
	} //while (step > 1)
}

void CDE_method::shell_sort_array_t(tFitness *array, int *list, int size)
{
	int done;
	int step, bound, i, j;
	int temp;

	for (i=0;i<size;i++)
	{
		list[i] = i;
	}

	step = size;  // array length
	while (step > 1) 
	{
		step /= 2;	//halve the step size
		do 
		{
			done   = 1;
			bound  = size - step;
			for (j = 0; j < bound; j++) 
			{
				i = j + step + 1;
				if (array[list[j]] < array[list[i-1]]) 	
				{
					temp      = list[i-1];
					list[i-1] = list[j];
					list[j]   = temp;
					done      = 0; // if a swap has been made we are not finished yet
				}  // if
			}  // for
		} while (done == 0);   // while
	} //while (step > 1)
}

double CDE_method::median(double *array, int size)
{
	double value = 0.0;

	int list[2000];
	for (int i=0;i<2000;i++)
	{
		list[i] = i;
	}
	shell_sort_array(array, list, size);
	
	if (size%2 == 1)
	{// the size of the array is odd.
		int media_index;
		media_index = size/2;
		value = array[list[media_index]];
	}
	else
	{// the size of the array is even.
		int a = (size-1)/2;
		int b = size/2;
		value = (array[list[a]]+array[list[b]])/2.0;
	}

	return value;
}

double CDE_method::mean_std(double *array, int size, double &stdv)
{
	double value = 0.0;

	assert(size > 0);

	for (int i=0;i<size;i++)
	{
		value += array[i];
	}
	value = value/((double)size);

	stdv = 0.0;
	for (int i =0;i<size;i++)
	{
		stdv += (array[i]-value)*(array[i]-value);
	}
	stdv = sqrt(stdv/((double)size-1));

	return value;
}

tFitness CDE_method::mean_std_t(tFitness *array, int size, tFitness &stdv)
{
	tFitness value = 0.0;

	assert(size > 0);

	for (int i=0;i<size;i++)
	{
		value += array[i];
	}
	value = value/((double)size);

	stdv = 0.0;
	for (int i =0;i<size;i++)
	{
		stdv += (array[i]-value)*(array[i]-value);
	}
	stdv = sqrt(stdv/((double)size-1));

	return value;
}

void CDE_method::display_result(long int gen)
{
	if(gen%10==0 || gen == 1)
	{
		//��ǰ��ø����������Ļ
		cout<<setw(5)<<gen;
		cout<<setw(8)<<evaluations;
		//��ʾԼ��������ֵ
		if (N_of_constr != 0)
		{			
			cout<<setw(15)<<best_individual.constr_violation;
		}
		cout.precision(10);					//�����������
		cout<<setw(20)<<MINIMIZE * parent_pop[best_index].fitness;
		cout<<setw(20)<<MINIMIZE * parent_pop[m_rnd.rndint(0, pop_size-1)].fitness;
		cout<<setw(20)<<MINIMIZE * parent_pop[worst_index].fitness<<endl;
	}
}

void CDE_method::report_result(long int gen, ofstream &file)
{
	if (gen < 5 )
	{
		//��ǰ��ø���������ļ�
		file<<setw(5)<<gen;
		file<<setw(10)<<evaluations;
		file.precision(20);			//�����������
		file<<setw(30)<<MINIMIZE * best_individual.fitness;
		file<<setw(30)<<best_individual.constr_violation;
		//file<<setw(30)<<best_individual.CR;
		file<<endl;
	}
	else
	{
		if ((gen % (max_iterations/25) == 0 ||
			gen >= max_iterations ||
			evaluations >= max_NFFEs) )
		{
			//��ǰ��ø���������ļ�
			file<<setw(5)<<gen;
			file<<setw(10)<<evaluations;
			file.precision(20);			//�����������
			file<<setw(30)<<MINIMIZE * best_individual.fitness;
			file<<setw(30)<<best_individual.constr_violation;
			//file<<setw(30)<<best_individual.CR;
			file<<endl;
		}
	}
}

void CDE_method::report_diversity(long int gen, ofstream &file)
{
	// calculate the diversity of the population
	double x_average[2000];
	for (int i=0;i<2000;i++)
	{
		x_average[i] = 0.0;
	}

	for (int i=0;i<N_of_x;i++)
	{
		for (int j=0;j<pop_size;j++)
		{
			x_average[i] += parent_pop[j].xreal[i];
		}
		x_average[i] = x_average[i]/((double)pop_size);
	}

	double diversity = 0.0;
	for (int i =0;i<pop_size;i++)
	{
		double temp = 0.0;
		double x;
		for (int j=0;j<N_of_x;j++)
		{
			x = parent_pop[i].xreal[j];
			temp += (x-x_average[j])*(x-x_average[j]);
		}
		diversity += sqrt(temp);
	}
	diversity = diversity/(((double)pop_size)*max_interval);

	pop_diversity = diversity;
	
	if (gen < 100 )
	{
		//��ǰ��ø���������ļ�
		file<<setw(5)<<gen;
		file<<setw(10)<<evaluations;
		file.precision(20);			//�����������
		file<<setw(30)<<diversity<<endl;
	}
	else
	{
		if ((gen % 25 == 0 ||
			gen >= max_iterations ||
			evaluations >= max_NFFEs) )
		{
			//��ǰ��ø���������ļ�
			file<<setw(5)<<gen;
			file<<setw(10)<<evaluations;
			file.precision(20);			//�����������
			file<<setw(30)<<diversity<<endl;
		}
	}
}

void CDE_method::report_parameter(long int gen, ofstream &file)
{	
	int counter = max_iterations/20;
	if (method_flag==7 || method_flag==8 )
	{
		if (gen <=10 || gen%counter == 0)
		{
			//��ǰ��ø���������ļ�
			file<<setw(5)<<gen;
			file<<setw(10)<<evaluations;
			file.precision(20);			//�����������
			file<<setw(30)<<"0";
			file<<setw(30)<<"0";
			file<<setw(30)<<"0";
			file<<endl;
		}
	}
	else
	{
		if (gen <=10 || gen%counter == 0)
		{
			//��ǰ��ø���������ļ�
			file<<setw(5)<<gen;
			file<<setw(10)<<evaluations;
			file.precision(20);			//�����������
			file<<setw(30)<<best_individual.CR;
			file<<setw(30)<<best_individual.F;
			file<<setw(30)<<"0";
			file<<endl;
		}
	}
}

void CDE_method::sort_pop(CIndividual *pop, int size)
{
	for (int i=0;i<size;i++)
	{
		sort_index[i] = i;
	}
	shell_sort_pop(pop, sort_index, size);
}

void CDE_method::Run_Optimizer(int run_no, double seed, int method, int strategy_ii)
{
	m_strategy = strategy_ii;

	if (m_strategy==1)
	{
		printf("Strategy:  rand/1/bin.\n");
	}
	
	switch(method) {
	case 1:
		printf("Algorithm: DE_method.\n");
		break;
	case 2:
		printf("Algorithm: rank_DE_method.\n");
		break;	
	default:
		printf("The method does not exist!\n");
		exit(0);
	}

	int i = 0; 
	
	method_flag = method;
	rnd_seed    = seed;

	int func_index;
	if (index_of_normal)	func_index = index_of_normal;

	ofstream SummaryFile;
	SummaryFile.open(".\\results\\summary.txt",ios::app);

	char f_name1[150];
	sprintf(f_name1,"%s%d%s",".\\results\\process\\process_",run_no+1,".txt");	
	ofstream file_process(f_name1);

	char f_name2[150];
	sprintf(f_name2,"%s%d%s",".\\results\\diversity\\diversity_",run_no+1,".txt");	
	ofstream file_diversity(f_name2);

	char f_name3[150];
	sprintf(f_name3,"%s%d%s",".\\results\\migration\\migration_",run_no+1,".txt");	
	ofstream file_migration(f_name3);

	clock_t start, finish;
	double time_consume;

	time_consume = 0.0;
	start = clock();						// starts the clock
	
	srand((unsigned)time(0));
	m_rnd.randomize(seed);

	evaluations = 0;						// reset the NFFEs

	init_variables();
	init_pop();
	evaluate_pop(parent_pop, pop_size);
	find_best_index(parent_pop, pop_size);
	best_individual = parent_pop[best_index];

	gen = 1;								// current generation number
	int counter=0;							// only for restarting the population

	report_result(gen, file_process);
	report_diversity(gen, file_diversity);
	report_parameter(gen, file_migration);

	int feasible_flag = 0;					// check to get the feasible individual first
	flag_precision = 0;						// check to arrive the required value
	/* -------- add different optimizer here -------- */
	while ( (evaluations < max_NFFEs) && (gen < max_iterations) )
	{
		finish = clock();					// time consuming of this generation
		time_consume = (double)(finish-start)/CLOCKS_PER_SEC;

		// report the results in the screen
		/* comment the following routines to save the time_consuming */
		//display_result(gen);
		
		// ADD YOUR OPTIMIZER HERE
		switch(method) {
		case 1:// DE method
			run_DE_method();
			break;
		case 2:// rank_DE method
			run_rank_DE_method();
			break;
		default:
			printf("The method selected does not exist.\n");
			exit(0);
		}

		gen++;

		find_best_index(parent_pop, pop_size);
		
		if (compare_ind(&parent_pop[best_index], &best_individual) != -1)
		{
			best_individual = parent_pop[best_index];
		}
		
		report_result(gen, file_process);
		report_diversity(gen, file_diversity);
		report_parameter(gen, file_migration);

		if (flag_precision == 0 && fabs(MINIMIZE*best_individual.fitness - known_optimal) < PRECISION && best_individual.feasible == 1)
		{
			flag_precision = 1;
			ofstream SummaryFile1;
			SummaryFile1.open(".\\results\\evaluations.txt",ios::app);				
			SummaryFile1<<setw(9)<<evaluations<<endl;
			SummaryFile1.close();
			//break;		// ������һ���ľ���ͳ����Ӧֵ���۴���ʱʹ��
		}
		if (N_of_constr != 0 && feasible_flag == 0 && best_individual.feasible == 1)
		{
			feasible_flag = 1;
			ofstream SummaryFile1;
			SummaryFile1.open(".\\results\\feasible_NFFEs.txt",ios::app);
			SummaryFile1<<setw(9)<<evaluations<<endl;
			SummaryFile1.close();
		}
	}
	/* -------- add different optimizer here -------- */
	
	printf("The total running time is %f s.\n", time_consume);
	if (index_of_normal != 0)
	{
		printf("The routine in optimizing f(%d) exits successfully.\n\n", index_of_normal);
	}

	SummaryFile.precision(15);
	SummaryFile<<setw(5)<<gen<<setw(9)<<evaluations
		<<setw(25)<<MINIMIZE * best_individual.fitness
		<<setw(8)<<time_consume
		<<setw(25)<<best_individual.constr_violation<<endl;
		
	file_process.close();
	file_diversity.close();
	file_migration.close();
	SummaryFile.close();
	
	return;
}

/************************************************************************/
/* For the original DE and jDE methods                                  */
/************************************************************************/
void CDE_method::run_DE_method()
{
	int    i, j;
	int    r1, r2, r3, r4, r5;
	double low, up;

	for (i=0;i<pop_size;i++)
	{
		// select five parents randomly
		do {
			r1 = m_rnd.rndint(0, pop_size-1);
		} while (r1 == i);
		do {
			r2 = m_rnd.rndint(0, pop_size-1);
		} while(r2 == i || r2 == r1);
		do { 
			r3 = m_rnd.rndint(0, pop_size-1);
		} while(r3 == i || r3 == r2 || r3 == r1);
		do { 
			r4 = m_rnd.rndint(0, pop_size-1);
		} while(r4 == i || r4 == r3 || r4 == r2 || r4 == r1);
		do { 
			r5 = m_rnd.rndint(0, pop_size-1);
		} while(r5 == i || r5 == r4 || r5 == r3 || r5 == r2 || r5 == r1);

		// parameter self-adaptation
		if (param_adaptation == 1)
		{
			child_pop[i].F  = parent_pop[i].F;
			child_pop[i].CR = parent_pop[i].CR;
			if (m_rnd.flip(0.1))
			{
				child_pop[i].F  = 0.1+m_rnd.randomperc()*0.9;
			}		
			if (m_rnd.flip(0.1))
			{
				child_pop[i].CR = m_rnd.rndreal(0.0, 1.0);
			}
			m_F  = child_pop[i].F;	// F
			m_CR = child_pop[i].CR;	// CR
		}
		else
		{
			m_CR = 0.9;
			m_F  = 0.5;//m_rnd.rndreal(0.1, 1.0);
		}

		/* mutation */
		for (j=0;j<N_of_x;j++)
		{
			low = lower_bounds[j];
			up  = upper_bounds[j];

			switch(m_strategy) {
			case 1:// DE/rand/1
				child_pop[i].xreal[j] = parent_pop[r1].xreal[j]
					+ m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]);
				break;
			case 2:// DE/rand/2
				child_pop[i].xreal[j]=parent_pop[r1].xreal[j]
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]))
					+(m_F*(parent_pop[r4].xreal[j]-parent_pop[r5].xreal[j]));
				break;
			case 3:// DE/current-to-best/1
				child_pop[i].xreal[j]=parent_pop[i].xreal[j]
					+(m_F*(parent_pop[best_index].xreal[j]-parent_pop[i].xreal[j]))
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]));
				break;
			case 4:// DE/current-to-best/2
				child_pop[i].xreal[j]=parent_pop[i].xreal[j]
					+(m_F*(parent_pop[best_index].xreal[j]-parent_pop[i].xreal[j]))
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]))
					+(m_F*(parent_pop[r4].xreal[j]-parent_pop[r5].xreal[j]));
				break;
			case 5:// DE/rand-to-best/1
				child_pop[i].xreal[j]=parent_pop[r1].xreal[j]
					+(m_F*(parent_pop[best_index].xreal[j]-parent_pop[r1].xreal[j]))
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]));
				break;
			case 6:// DE/rand-to-best/2
				child_pop[i].xreal[j]=parent_pop[r1].xreal[j]
					+(m_F*(parent_pop[best_index].xreal[j]-parent_pop[r1].xreal[j]))
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]))
					+(m_F*(parent_pop[r4].xreal[j]-parent_pop[r5].xreal[j]));
				break;
			default:
				printf("The scheme in DE method does not exist.\n");
				exit(0);
			}

			if (child_pop[i].xreal[j] < low || child_pop[i].xreal[j] > up)
			{
				child_pop[i].xreal[j] = m_rnd.rndreal(low, up);
			}
		}
	
		/* crossover */
		int j_rand = m_rnd.rndint(0, N_of_x-1);
		for (j=0;j<N_of_x;j++)
		{
			if (1)//m_strategy != 6)
			{
				if (m_rnd.rndreal(0,1)<m_CR || j==j_rand)
				{
					child_pop[i].xreal[j] = child_pop[i].xreal[j];
				}
				else
				{
					child_pop[i].xreal[j] = parent_pop[i].xreal[j];
				}
			}
			else
			{// only for DE/current-to-rand/1 strategy (crossover is not used)
				child_pop[i].xreal[j] = child_pop[i].xreal[j];
			}
		}
	}

	evaluate_pop(child_pop, pop_size);

	for (i=0;i<pop_size;i++)
	{
		int flag = compare_ind(&child_pop[i], &parent_pop[i]);
		if (flag != -1)
		{
			parent_pop[i] = child_pop[i];
		}
	}
}

/************************************************************************/
/* DE with probability-based mutation operators                         */
/************************************************************************/
/* Calculate the probabilities of each individual based on its fitness  */
void CDE_method::normalize_reward(tFitness * S_reward, int S_size)
{
	int i;
	/* normalize the reward */
	tFitness max_rr = S_reward[0];
	for (i=1;i<S_size;i++)
	{
		if (S_reward[i] > max_rr)
		{
			max_rr = S_reward[i];
		}
	}
	if (max_rr != 0.0)
	{
		for (i=0;i<S_size;i++)
		{
			S_reward[i] = S_reward[i]/max_rr;
		}
	}
}

void CDE_method::calculate_weights(tFitness *S_reward, int S_size, double *weights, int flag)
{
	int i;

	/* sort the rewards in descending order */
	int weight_index[1000];
	shell_sort_array_t(S_reward, weight_index, S_size);

	if (1==flag)
	{// rank-based weights: linear model
		for (i=0;i<S_size;i++)
		{
			weights[weight_index[i]] = (double)(i+0)/S_size;					// linear migration model
		}
	}
	else if (2==flag)
	{// rank-based weights: quadratic model
		for (i=0;i<S_size;i++)
		{
			weights[weight_index[i]] = pow( (double)(i+0)/S_size, 2.0 );		// quadratic migration model
		}
	}
	else if (3==flag)
	{// rank-based weights: cosine model
		for (i=0;i<S_size;i++)
		{
			weights[weight_index[i]] = 0.5*(-cos((i+0)*PI/S_size) + 1.0);		// sinusoidal migration model
		}
	}
	else
	{
		printf("The weights calculating is wrong.\n");
		exit(0);
	}
}

void CDE_method::rank_vector_selection(int base_index, int *index, int size, const double *probability)
{
	assert(size>=2 && size<=5);

	int i = base_index;

	if (2==size)
	{// select two different vectors
		do {
			index[0] = m_rnd.rndint(0, pop_size-1);
		} while (index[0]==i || !m_rnd.flip(probability[index[0]]));
		do {
			index[1] = m_rnd.rndint(0, pop_size-1);
		} while(index[1]==i || index[1]==index[0]);
	}
	else if (3==size)
	{// select three different vectors
		do {
			index[0] = m_rnd.rndint(0, pop_size-1);
		} while (index[0]==i || !m_rnd.flip(probability[index[0]]));
		do {
			index[1] = m_rnd.rndint(0, pop_size-1);
		} while(index[1]==i || index[1]==index[0] || !m_rnd.flip(probability[index[1]]));		
		do {
			index[2] = m_rnd.rndint(0, pop_size-1);
		} while(index[2]==i || index[2]==index[1] || index[2]==index[0] || m_rnd.flip(probability[index[2]]));
	}
	else if (4==size)
	{// select five different vectors
		do {
			index[0] = m_rnd.rndint(0, pop_size-1);
		} while (index[0]==i || !m_rnd.flip(probability[index[0]]));
		do {
			index[1] = m_rnd.rndint(0, pop_size-1);
		} while(index[1]==i || index[1]==index[0] || !m_rnd.flip(probability[index[1]]));		
		do {
			index[2] = m_rnd.rndint(0, pop_size-1);
		} while(index[2]==i || index[2]==index[1] || index[2]==index[0]);			
		do {
			index[3] = m_rnd.rndint(0, pop_size-1);
		} while(index[3]==i || index[3]==index[2] || index[3]==index[1] || index[3]==index[0] || !m_rnd.flip(probability[index[3]]));
	}
	
	else// if (5==size)
	{// select five different vectors
		do {
			index[0] = m_rnd.rndint(0, pop_size-1);
		} while (index[0]==i || !m_rnd.flip(probability[index[0]]));
		do {
			index[1] = m_rnd.rndint(0, pop_size-1);
		} while(index[1]==i || index[1]==index[0] || !m_rnd.flip(probability[index[1]]));		
		do {
			index[2] = m_rnd.rndint(0, pop_size-1);
		} while(index[2]==i || index[2]==index[1] || index[2]==index[0]);			
		do {
			index[3] = m_rnd.rndint(0, pop_size-1);
		} while(index[3]==i || index[3]==index[2] || index[3]==index[1] || index[3]==index[0] || !m_rnd.flip(probability[index[3]]));
		do {
			index[4] = m_rnd.rndint(0, pop_size-1);
		} while(index[4]==i || index[4]==index[3] || index[4]==index[2] || index[4]==index[1] || index[4]==index[0]);
	}
}

void CDE_method::run_rank_DE_method()
{
	int    i, j;
	int    r1, r2, r3, r4, r5;
	double low, up;

	/************************************************************************/
	/* Calculate the probabilities of each individual based on its fitness  */
	/************************************************************************/
	double   probability[population_size];
	tFitness fitness[population_size];
	for (i=0; i<pop_size; i++)
	{
		fitness[i] = parent_pop[i].fitness;
	}
	calculate_weights(fitness, pop_size, probability, weights_flag);
	int indices[10] = {0};
	/************************************************************************/

	for (i=0;i<pop_size;i++)
	{
		// select five parents randomly
		rank_vector_selection(i, indices, 3, probability);////////////////////////////////
		r1 = indices[0];	r2 = indices[1];	r3 = indices[2];	
		r4 = indices[3];	r5 = indices[4];

		// parameter self-adaptation
		if (param_adaptation == 1)
		{
			child_pop[i].F  = parent_pop[i].F;
			child_pop[i].CR = parent_pop[i].CR;
			if (m_rnd.flip(0.1))
			{
				child_pop[i].F  = 0.1+m_rnd.randomperc()*0.9;
			}		
			if (m_rnd.flip(0.1))
			{
				child_pop[i].CR = m_rnd.rndreal(0.0, 1.0);
			}
			m_F  = child_pop[i].F;	// F
			m_CR = child_pop[i].CR;	// CR
		}
		else
		{
			m_CR = 0.9;
			m_F  = 0.5;//m_rnd.rndreal(0.1, 1.0);
		}

		/* mutation */
		for (j=0;j<N_of_x;j++)
		{
			low = lower_bounds[j];
			up  = upper_bounds[j];

			switch(m_strategy) {
			case 1:// DE/rand/1
				child_pop[i].xreal[j] = parent_pop[r1].xreal[j]
					+ m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]);
				break;
			case 2:// DE/rand/2
				child_pop[i].xreal[j]=parent_pop[r1].xreal[j]
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]))
					+(m_F*(parent_pop[r4].xreal[j]-parent_pop[r5].xreal[j]));
				break;
			case 3:// DE/current-to-best/1
				child_pop[i].xreal[j]=parent_pop[i].xreal[j]
					+(m_F*(parent_pop[best_index].xreal[j]-parent_pop[i].xreal[j]))
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]));
				break;
			case 4:// DE/current-to-best/2
				child_pop[i].xreal[j]=parent_pop[i].xreal[j]
					+(m_F*(parent_pop[best_index].xreal[j]-parent_pop[i].xreal[j]))
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]))
					+(m_F*(parent_pop[r4].xreal[j]-parent_pop[r5].xreal[j]));
				break;
			case 5:// DE/rand-to-best/1
				child_pop[i].xreal[j]=parent_pop[r1].xreal[j]
					+(m_F*(parent_pop[best_index].xreal[j]-parent_pop[r1].xreal[j]))
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]));
				break;
			case 6:// DE/rand-to-best/2
				child_pop[i].xreal[j]=parent_pop[r1].xreal[j]
					+(m_F*(parent_pop[best_index].xreal[j]-parent_pop[r1].xreal[j]))
					+(m_F*(parent_pop[r2].xreal[j]-parent_pop[r3].xreal[j]))
					+(m_F*(parent_pop[r4].xreal[j]-parent_pop[r5].xreal[j]));
				break;
			default:
				printf("The scheme in DE method does not exist.\n");
				exit(0);
			}

			if (child_pop[i].xreal[j] < low || child_pop[i].xreal[j] > up)
			{
				child_pop[i].xreal[j] = m_rnd.rndreal(low, up);
			}
		}
	
		/* crossover */
		int j_rand = m_rnd.rndint(0, N_of_x-1);
		for (j=0;j<N_of_x;j++)
		{
			if (1)//m_strategy != 6)
			{
				if (m_rnd.rndreal(0,1)<m_CR || j==j_rand)// || m_rnd.flip(1.0-probability[i]))
				{
					child_pop[i].xreal[j] = child_pop[i].xreal[j];
				}
				else
				{
					child_pop[i].xreal[j] = parent_pop[i].xreal[j];
				}
			}
			else
			{// only for DE/current-to-rand/1 strategy (crossover is not used)
				child_pop[i].xreal[j] = child_pop[i].xreal[j];
			}
		}
	}

	evaluate_pop(child_pop, pop_size);

	for (i=0;i<pop_size;i++)
	{
		int flag = compare_ind(&child_pop[i], &parent_pop[i]);
		if (flag != -1)
		{
			parent_pop[i] = child_pop[i];
		}
	}
}

