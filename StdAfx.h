// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__DD406DC6_1ADE_41F7_AD36_39F133CED382__INCLUDED_)
#define AFX_STDAFX_H__DD406DC6_1ADE_41F7_AD36_39F133CED382__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


// TODO: reference additional headers your program requires here

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

// --------------------include files--------------------
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <windows.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <assert.h>


// --------------------Global variables--------------------
#define E	exp(1.0)
#define PI	acos(-1.0)
#define INF	1.0e99
#define EPS 1.0e-14		// the precision of the constraints

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

#define DBL_MAX 1.7976931348623158e+308 

typedef long double tFitness;


#endif // !defined(AFX_STDAFX_H__DD406DC6_1ADE_41F7_AD36_39F133CED382__INCLUDED_)
