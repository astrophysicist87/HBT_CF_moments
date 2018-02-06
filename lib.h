#ifndef LIB_H
#define LIB_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_integration.h>

using namespace std;

const complex<double> i(0, 1);

inline complex<double> cot(complex<double> x){return (cos(x)/sin(x));}
inline complex<double> csc(complex<double> x){return (1.0/sin(x));}
inline complex<double> sec(complex<double> x){return (1.0/cos(x));}

void linspace(vector<double> & x, double a, double b);

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

#endif
