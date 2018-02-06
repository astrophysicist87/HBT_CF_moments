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
#include <cfloat>

//#include "declarations.h"
#include "gauss.h"
//#include "qng_1d_vec.h"

using namespace std;

	const double PI = 3.14159265358979323846264338327950;
	const int order = 43;
	const int half_order = (order-1)/2;
	const double x_lower=-6.25, x_upper=6.25;
	const double x_interval = (x_upper-x_lower)/double(order);
	const double half_length = (x_upper-x_lower)*0.5;
	const double center = (x_upper+x_lower)*0.5;
	const double overall_scale_factor = 10.;
	double wts[order];
	double xpts[order];
	double wts_mod[order];
	double xpts_mod[order];

void test_function(double x, double * result);

int
main ()
{
	time_t start, end;

	//cout << "The largest double is " << std::numeric_limits<double>::infinity() << endl;

	cout << "The largest double is " << DBL_MAX << endl;

	time(&start);

	for (int i = 0;  i <= 20; i++) {
		double my_answer;
		test_function(double(i)/20., &my_answer);
		cout << "i = " << i << " \t f = " << my_answer << endl;
	}

	time(&end);
	cout << "The whole process took " << difftime(end,start) << " seconds" << endl;

	return (0);
}

void test_function(double x, double * result)
{
	*result = exp(-x*x+x-1.)*cos(x);
}

//End of file
