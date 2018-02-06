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
	const int order = 43;
	const int half_order = (order-1)/2;
	const double x_lower=-6.5, x_upper=6.5;
	const double x_interval = (x_upper-x_lower)/double(order);
	const double half_length = (x_upper-x_lower)*0.5;
	const double center = (x_upper+x_lower)*0.5;
	const double overall_scale_factor = 10.;
	double wts[order];
	double xpts[order];
	double wts_mod[order];
	double xpts_mod[order];

double test_function(double x);
double qng_1d_vec (vector<double>* vector_ptr, double a1, double b1);
void generate_vector(double (*function_original) (double), vector<double>* vector_ptr, double a1, double b1);
double qng_1d_vec_new (vector<double>* vector_ptr);
void generate_vector_new(double (*function_original) (double), vector<double>* vector_ptr);
void new_points_and_weights(int mod_index);

int
main ()
{
	time_t start, end;

	time(&start);
	gauss(order, 3, -1, 1, xpts, wts);
	new_points_and_weights(2);	//1: semi-infinite interval
					//2: infinite interval

/*	for (int i = 0; i < order; i++) {
		cout << "xpts[" << i << "] = " << xpts[i] << " and wts[" << i << "] = " << wts[i] << endl;
	}
*/
	vector<double>* vec_ptr;
	vec_ptr = new vector<double> (order);
	vector <double>* new_vec_ptr;
	new_vec_ptr = new vector<double> (order);

	generate_vector(&test_function, vec_ptr, x_lower, x_upper);
	generate_vector_new(&test_function, new_vec_ptr);

	double result = qng_1d_vec(vec_ptr, x_lower, x_upper);
	double result_new = qng_1d_vec_new(new_vec_ptr);

	cout << setprecision(20) << "The integral = " << result << endl;
	cout << setprecision(20) << "The new integral = " << result_new << endl;
	cout << setprecision(20) << "The exact integral = " << sqrt(PI) << endl;

	time(&end);
	cout << "The whole process took " << difftime(end,start) << " seconds" << endl;

	return (0);
}

double test_function(double x)
{
	const double result = exp(-x*x);

	return (result);
}

double
qng_1d_vec (vector<double>* vector_ptr, double a1, double b1)  	//lower and upper limits for x integration
{
	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	double result_int = 0.;

	for (int i = 0; i <= half_order; i++) {
		const double fval = (*vector_ptr).at(order-1-i) + (*vector_ptr).at(i);
			//this part sums up all function values at different abscissae which receive the same weights
		double factor = 1.;
		if (i == half_order) {factor /= 2.;}
		result_int += fval * wts[i] * factor;
	}
	result_int *= half_length1;

	return (result_int);
}

void generate_vector(double (*function_original) (double), vector<double>* vector_ptr, double a1, double b1)
{
	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	double result_int = 0.;

	for (int i = 0; i <= order-1; i++) {
		double abscissa1;
		//if (i < half_order) {abscissa1 = -half_length1 * xpts[i];}
		//else if (i == half_order) {abscissa1 = 0.;}
		//else {abscissa1 = half_length1 * xpts[order-1-i];}
		abscissa1 = half_length1 * xpts[i];
		const double fval = function_original(center1 + abscissa1);
		(*vector_ptr).at(i) = fval;
	}
}

double
qng_1d_vec_new (vector<double>* vector_ptr)
{
	const double midpoint = 1.;
	double result_int = 0.;

	for (int i = 0; i <= order-1; i++) {
		result_int += wts_mod[i] * (*vector_ptr).at(i);
	}

	return (result_int);
}

void generate_vector_new(double (*function_original) (double), vector<double>* vector_ptr)
{
	for (int i = 0; i <= order-1; i++) {
		(*vector_ptr).at(i) = function_original(xpts_mod[i]);
	}
}

void new_points_and_weights(int mod_index)
{
	const double midpoint = 1.;
	if (mod_index == 1) {  //half-infinite
		for (int i = 0; i <= order - 1; i++) {
			xpts_mod[i] = midpoint * (1. + xpts[i]) / (1. - xpts[i]);
			wts_mod[i] = 2. * midpoint * wts[i] / ( (1. - xpts[i]) * (1. - xpts[i]) );
		}
	}
	else {  //infinite
		for (int i = 0; i <= order - 1; i++) {
			xpts_mod[i] = midpoint * xpts[i] / (1. - xpts[i] * xpts[i]);
			wts_mod[i] = midpoint * wts[i] * (1. + xpts[i] * xpts[i]) / ( (1. - xpts[i] * xpts[i]) * (1. - xpts[i] * xpts[i]) );
		}
	}
}

//End of file
