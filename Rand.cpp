// Rand.cpp: implementation of the CRand class.
//
//////////////////////////////////////////////////////////////////////

// #include "StdAfx.h"
#include <math.h>
#include <assert.h>
#include "Rand.h"
const double PI = atan(1.0)*4;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

int    CRand::jrand = 0;
double CRand::oldrand[55] = {0};
double CRand::rndx1 = 0;
double CRand::rndx2 = 0;
int    CRand::rndcalcflag = 1;

CRand::CRand()
{
}

CRand::~CRand()
{
}

// Get seed number for random and start it up
void CRand::randomize(double seed)
{
	int j1;
	for(j1=0; j1<=54; j1++)
    {
        oldrand[j1] = 0.0;
    }
	jrand=0;
	warmup_random (seed);
    return;
}

// Get randomize off and running
void CRand::warmup_random (double seed)
{
    int j1, ii;
    double new_random, prev_random;
    oldrand[54] = seed;
    new_random = 0.000000001;
    prev_random = seed;
    for(j1=1; j1<=54; j1++)
    {
        ii = (21*j1)%54;
        oldrand[ii] = new_random;
        new_random = prev_random-new_random;
        if(new_random<0.0)
        {
            new_random += 1.0;
        }
        prev_random = oldrand[ii];
    }
    advance_random ();
    advance_random ();
    advance_random ();
    jrand = 0;
    return;
}

// Create next batch of 55 random numbers
void CRand::advance_random ()
{
    int j1;
    double new_random;
    for(j1=0; j1<24; j1++)
    {
        new_random = oldrand[j1]-oldrand[j1+31];
        if(new_random<0.0)
        {
            new_random = new_random+1.0;
        }
        oldrand[j1] = new_random;
    }
    for(j1=24; j1<55; j1++)
    {
        new_random = oldrand[j1]-oldrand[j1-24];
        if(new_random<0.0)
        {
            new_random = new_random+1.0;
        }
        oldrand[j1] = new_random;
    }
}

// Fetch a single random number between 0.0 and 1.0
double CRand::randomperc()
{
	jrand++;
	if(jrand>=55)
    {
        jrand = 1;
        advance_random();
    }
    return((double)oldrand[jrand]);
}

// Fetch a single random integer between low and high including the bounds
int CRand::rndint (int low, int high)
{
    int res;
    if (low >= high)
    {
        res = low;
    }
    else
    {
        res = low + (int)(randomperc()*(high-low+1));
        if (res > high)
        {
            res = high;
        }
    }
    return (res);
}

// Fetch a single random real number between low and high including the bounds
double CRand::rndreal (double low, double high)
{
    return (low + (high-low)*randomperc());
}

void CRand::initrandomnormaldeviate()
{
	rndcalcflag = 1;
}

double CRand::randomnormaldeviate()
{
	double t; 
 
    if(rndcalcflag) 
    { 
        rndx1 = sqrt(- 2.0*log((double) randomperc())); 
        t = 6.2831853072 * (double) randomperc(); 
        rndx2 = sin(t); 
        rndcalcflag = 0; 
        return(rndx1 * cos(t)); 
    } 
    else 
    { 
        rndcalcflag = 1; 
        return(rndx1 * rndx2); 
    }
}

double CRand::noise(double mu,double sigma)
{
	return((randomnormaldeviate()*sigma) + mu); 
}

double CRand::cauchy(double median, double formFactor)
{
	assert( formFactor > 0. ); 
	double u, v; 
	do { 
			double U1 = (double) randomperc();
			double U2 = (double) randomperc();
			u = 2.0 * U1 - 1.0; 
			v = 2.0 * U2 - 1.0; 
		}while ((u*u+v*v>1.0) || (u==0.0&&v==0.0)); 

	if (u!=0) 
		return (median + formFactor * (v/u)); 
	else 
		return (median);
}

double CRand::gaussian(double median, double stdev)
{
	assert( stdev > 0. ); 
	double u,v,x; 
	do { 
			double U1 = (double) randomperc();
			double U2 = (double) randomperc();
			u = 2.0 * U1 - 1.0; 
			v = 2.0 * U2 - 1.0; 
			x = u * u + v * v; 
		}while ( x >= 1.0 || x == 0 ); 
	return median + stdev * u * sqrt( -2. * log( x ) / x );
}

void CRand::bivariate_gaussian(double mean_x, double mean_y, double sigma_x, 
							   double sigma_y, double rho, double *x, double *y)
{
	double u, v, r2, scale;

	do
	{
		/* choose x,y in uniform square (-1,-1) to (+1,+1) */
		u = -1 + 2 * (double) randomperc();
		v = -1 + 2 * (double) randomperc();

		/* see if it is in the unit circle */
		r2 = u * u + v * v;
    } while (r2 > 1.0 || r2 == 0);

	scale = sqrt (-2.0 * log (r2) / r2);

	*x = mean_x + sigma_x * (rho * u + sqrt(1 - rho*rho) * v) * scale;
	*y = mean_y + sigma_y * u * scale;

//	if (u!=0) 
//		*y = mean_y + sigma_y * (v/u);
//	else 
//		*y = mean_y;
}

double CRand::log_normal(double median, double stdv)
{
	assert( stdv > 0. && median > 0.); 

	double log_stdv;
	double log_mean;

	double rst;

	log_stdv = sqrt(log(1+(stdv/median)*(stdv/median)));
	log_mean = log(median)-0.5*(1+(stdv/median)*(stdv/median));

	double U1 = (double) randomperc();
	double U2 = (double) randomperc();

	rst = sqrt(-2*log(U1))*sin(2*PI*U2);
	rst = log_mean + rst*log_stdv;
	rst = exp(rst);

	return rst;
}

int CRand::flip(double p)
{
	if (randomperc() <= p)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}