//Program: multidim_integ_test
//Author: Christopher Plumberg
//Date: November 26, 2012
//Modified from integration code written by Dick Furnstahl
//
//Comments: this algorithm only uses constant integration limits!
//May need to generalize to variable boundaries in future.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>

using namespace std;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

struct phys_params
{
	double alpha, beta;
};

//function prototype(s)
phys_params get_parameters();
double S_function (double x, double y, double eta, void *);
double eta_t(double x, double y, void *);
double x_mean(double x, double y, double eta, void *);
double y_mean(double x, double y, double eta, void *);
double eps_3_centershift_real(double x, double y, double eta, void *);
double eps_3_centershift_imag(double x, double y, double eta, void *);

/************************************************************************/
//Note that inclusion of header comes AFTER definition of phys_params struct!!!
/************************************************************************/
#include "gsl_3d_int.h"
/************************************************************************/

const double PI = 3.1415926535897932384626433;

int
main ()
{
  //Define upper and lower limits of integration
  double limit=10.;
  double xlo=-limit;
  double xup=limit;
  double ylo=-limit;
  double yup=limit;
  double etalo=-limit;
  double etaup=limit;
  double norm;
  phys_params initial_params = get_parameters();
  double plumbergtest;
  double eps_3_real;
  double eps_3_imag;
  int index1 = 10;  //number of steps of v_3_bar
  int index2 = 4;   //number of steps of eps_3_bar
  double array[4];
  double arrayre[4];
  double arrayim[4];

  gsl_set_error_handler_off();

  norm = gsl_3d_int_driver(&S_function, &initial_params, xlo, xup, ylo, yup, etalo, etaup);

  cout << "The answer is " << setprecision(25) << norm << endl;

  return (0);
}

/************************************************************************/

phys_params get_parameters()
{
	struct phys_params inits;

	//Define parameters for function to be integrated here
	inits.alpha=1.;
	inits.beta=1.;

	return (inits);
}

/************************************************************************/

double S_function (double x, double y, double eta, void *params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	double alpha = params.alpha;
	double beta = params.beta;

	return (exp(-alpha*(x*x+y*y+eta*eta)+beta));
}

//End of file
