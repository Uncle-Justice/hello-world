// Rand.h: interface for the CRand class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RAND_H__0003D1ED_F550_44C8_85BE_39361D1A3BAE__INCLUDED_)
#define AFX_RAND_H__0003D1ED_F550_44C8_85BE_39361D1A3BAE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CRand  
{
public:
	CRand();
	virtual ~CRand();

public:
	// Variable declarations for the random number generator
	static double oldrand[55];						// Array of 55 random numbers
	static int jrand;								// current random number 
	static double rndx1, rndx2;						// used with random normal deviate
	static int rndcalcflag;							// used with random normal deviate 

public:
	// Function declarations for the random number generator
	void   randomize(double seed);
	void   warmup_random (double seed);
	void   advance_random (void);
	double randomperc(void);
	int    rndint (int low, int high);				// integer value 0~n-1
	double rndreal (double low, double high);		// real value
	void   initrandomnormaldeviate();				// initialization routine for randomnormaldeviate
	double noise(double mu, double sigma);			// normal noise with specified mean & std dev: mu & sigma
	double randomnormaldeviate();					// random normal deviate after ACM algorithm 267 / Box-Muller Method
	
	double gaussian(double median, double stdv);	// Gaussian random number
	double log_normal(double median, double stdv);	// Log-normal random number
	double cauchy(double median, double factor);	// Cauchy random number
	int    flip(double);

	/* 
		This function computes the probability density p(x,y) at (x,y) for a bivariate gaussian distribution 
		with standard deviations sigma_x, sigma_y and correlation coefficient rho, using the formula given above.  
	*/
	void   bivariate_gaussian(double mean_x, double mean_y,
		                      double sigma_x, double sigma_y, 
		                      double rho,
                              double *x, double *y);// Bivariate Gaussian random number
};

#endif // !defined(AFX_RAND_H__0003D1ED_F550_44C8_85BE_39361D1A3BAE__INCLUDED_)
