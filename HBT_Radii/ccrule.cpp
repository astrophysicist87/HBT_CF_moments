//Program: R2_s_calc
//Author: Christopher Plumberg
//Date: November 26, 2012
//Modified from integration code written by Dick Furnstahl
//
//Comments: this algorithm only uses constant integration limits!
//May need to generalize to variable boundaries in future.
//
//Does not implement tau integration, since it is trivially (done analytically) in current model

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>
#include <complex>
#include <time.h>
#include <vector>

//#include "declarations.h"
#include "gauss.h"
//#include "qng_1d_vec.h"

using namespace std;

	const double PI = 3.14159265358979323846264338327950;
	const int order = 87;
	const int half_order = (order-1)/2;
	const double x_lower=-PI, x_upper=PI;
	const double x_interval = (x_upper-x_lower)/double(order);
	const double half_length = (x_upper-x_lower)*0.5;
	const double center = (x_upper+x_lower)*0.5;
	double wts[order];
	double abscissae[order];

double test_function(double x);

int
main ()
{
	time_t start, end;

	time(&start);
	gauss(order, 3, -1, 1, abscissae, wts);

	time(&end);
	cout << "The whole process took " << difftime(end,start) << " seconds" << endl;

	return (0);
}

double test_function(double x)
{
	//i want to find the sine and cosine transforms of the following function:
	const double result = exp(-cos(x)*cos(x));

	return (result);
}

//End of file
