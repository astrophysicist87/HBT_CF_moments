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
	double M_perp, T_0, eta_0, Y_rapidity;
	double Rad, eta_f, Phi_k, tau_f, K_perp;
	double Delta_tau, Delta_eta;

	double v_3_bar, eps_3_bar, psi_3_bar;

	double xavg, yavg;
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

  // open the output file stream
  ofstream output ("eps3_psi3bar_120.dat");	// save data
  ofstream outputre ("eps3_psi3bar_120_re.dat");
  ofstream outputim ("eps3_psi3bar_120_im.dat");

  output << "# v_3_bar \t eps_3_bar = 0 \t eps_3_bar = 1/8 \t eps_3_bar = 1/4 \t eps_3_bar = 3/8" << endl;
  outputre << "# v_3_bar \t eps_3_bar = 0 \t eps_3_bar = 1/8 \t eps_3_bar = 1/4 \t eps_3_bar = 3/8" << endl;
  outputim << "# v_3_bar \t eps_3_bar = 0 \t eps_3_bar = 1/8 \t eps_3_bar = 1/4 \t eps_3_bar = 3/8" << endl;

/*  cout << S_function(-1.,-1.,-1.,&initial_params) << endl;
  cout << eta_t(-1.,-1.,&initial_params) << endl;
  cout << x_mean(-1.,-1.,-1.,&initial_params) << endl;
  cout << y_mean(-1.,-1.,-1.,&initial_params) << endl;
  cout << eps_3_centershift_real(-1.,-1.,-1.,&initial_params) << endl;
  cout << eps_3_centershift_imag(-1.,-1.,-1.,&initial_params) << endl;*/
//  cin >> plumbergtest;

  for (int i = 1; i <= index1; i++)
  {
	initial_params.v_3_bar = double(i-1)/(2.*double(index1));
	for (int j = 1; j <= index2; j++)
	{
		initial_params.eps_3_bar = double(j-1)/(2.*double(index2));
		//norm = gsl_3d_int_driver(&goofy_extra_function, initial_params, xlo,
//			xup, ylo, yup, etalo, etaup);
		norm = gsl_3d_int_driver(&S_function, &initial_params, xlo,
					xup, ylo, yup, etalo, etaup);
		initial_params.xavg = (1./norm) * gsl_3d_int_driver(&x_mean, &initial_params, xlo,
					xup, ylo, yup, etalo, etaup);
		initial_params.yavg = (1./norm) * gsl_3d_int_driver(&y_mean, &initial_params, xlo,
					xup, ylo, yup, etalo, etaup);
		eps_3_real = (1./norm) * gsl_3d_int_driver(&eps_3_centershift_real, &initial_params, xlo,
					xup, ylo, yup, etalo, etaup);
		eps_3_imag = (1./norm) * gsl_3d_int_driver(&eps_3_centershift_imag, &initial_params, xlo,
					xup, ylo, yup, etalo, etaup);
		array[j-1] = sqrt(eps_3_real*eps_3_real + eps_3_imag*eps_3_imag);
		arrayre[j-1] = eps_3_real;
		arrayim[j-1] = eps_3_imag;
		cout << setprecision(8)  << "with v_3_bar =  " << initial_params.v_3_bar
		     << " and eps_3_bar =  " << initial_params.eps_3_bar << endl
		     << "we have norm = " << norm << ", " << "eps_3_real = " << eps_3_real
		     << " and eps_3_imag = " << eps_3_imag << endl
		     << "xavg = " << initial_params.xavg << " and yavg = " << initial_params.yavg << endl;
	}
	output << initial_params.v_3_bar << "\t" << array[0] << "\t" << array[1] << "\t" << array[2] << "\t" << array[3] << endl;
	outputre << initial_params.v_3_bar << "\t" << arrayre[0] << "\t" << arrayre[1] << "\t" << arrayre[2] << "\t" << arrayre[3] << endl;
	outputim << initial_params.v_3_bar << "\t" << arrayim[0] << "\t" << arrayim[1] << "\t" << arrayim[2] << "\t" << arrayim[3] << endl;
  }

//  cout << setprecision(25) << scientific << "norm is " << norm << endl
//       << "xavg is " << xavg << endl
//       << "yavg is " << yavg << endl;

  cout << "Data stored in eps3_psi3bar_120.dat, eps3_psi3bar_120_real.dat and eps3_psi3bar_120_imag.dat" << endl;
  output.close ();
  outputre.close ();
  outputim.close ();

  return (0);
}

/************************************************************************/

phys_params get_parameters()
{
	struct phys_params inits;

	//Define parameters for function to be integrated here
	inits.M_perp = 1.;
	inits.T_0 = 1.;
	inits.eta_0 = 0.;
	inits.Y_rapidity = 0.;
	inits.Rad = 1.;
	inits.eta_f = 1.;
	inits.Phi_k = 0.;
	inits.tau_f = 0.;
	inits.K_perp = 1.;
	inits.eta_0 = 0.;
	inits.Delta_tau = 1.;
	inits.Delta_eta = 1.;

	inits.v_3_bar = 0.;
	inits.eps_3_bar = 0.;
	inits.psi_3_bar = 8.*PI/12.;

	inits.xavg = 0.;
	inits.yavg = 0.;

	return (inits);
}

/************************************************************************/

double S_function (double x, double y, double eta, void *params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	double M_perp = params.M_perp;
	double T_0 = params.T_0;
	double eta_0 = params.eta_0;
	double Y_rapidity = params.Y_rapidity;
	double Rad = params.Rad;
	double Phi_k = params.Phi_k;
	double tau_f = params.tau_f;
	double K_perp = params.K_perp;
	double Delta_tau = params.Delta_tau;
	double Delta_eta = params.Delta_eta;

	double v_3_bar = params.v_3_bar;
	double eps_3_bar = params.eps_3_bar;
	double psi_3_bar = params.psi_3_bar;

	//Need to express everything in Cartesian coordinates
	//for greatest simplicity in integration
	double r = sqrt(x*x+y*y);
	double cos_phi = x/r;
	double sin_phi = y/r;
	double cos_3_phi = 4.*cos_phi*cos_phi*cos_phi - 3.*cos_phi;
	double sin_3_phi = 3.*sin_phi - 4.*sin_phi*sin_phi*sin_phi;

	//Expression for function is long and complicated,
	//so break it up into constituent terms
	double term1, term2, term3, term4;
	term1 = (eta-eta_0)*(eta-eta_0) / (2.*Delta_eta*Delta_eta);
	//term1 = 0.;
	term2 = (M_perp/T_0) * cosh(eta - Y_rapidity) * cosh(eta_t(x, y, &params));
	//term2 = 0.;
	term3 = (K_perp/T_0) * (cos_phi * cos(Phi_k) + sin_phi * sin(Phi_k)) * sinh(eta_t(x, y, &params));
	//term3 = 0.;
	term4 = (r*r)/(2.*Rad*Rad) * (1.-2.*eps_3_bar*( cos_3_phi * cos(3.*psi_3_bar) + sin_3_phi * sin(3.*psi_3_bar) ));
	//term4 = 0.;

	return (exp(-term1 - term2 + term3 - term4));
	//return (exp(-(x*x+y*y+eta*eta)+1.));
}

double eta_t(double x, double y, void *params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	double eta_f = params.eta_f;
	double Rad = params.Rad;

	double v_3_bar = params.v_3_bar;
	double eps_3_bar = params.eps_3_bar;
	double psi_3_bar = params.psi_3_bar;

	//Need to express everything in Cartesian coordinates
	//for greatest simplicity in integration
	double r = sqrt(x*x+y*y);
	double cos_phi = x/r;
	double sin_phi = y/r;
	double cos_3_phi = 4.*cos_phi*cos_phi*cos_phi - 3.*cos_phi;
	double sin_3_phi = 3.*sin_phi - 4.*sin_phi*sin_phi*sin_phi;

	return (eta_f * (r/Rad) * (1.+2.*v_3_bar*( cos_3_phi * cos(3.*psi_3_bar) + sin_3_phi * sin(3.*psi_3_bar) )));
}

double x_mean(double x, double y, double eta, void *params_ptr)
{
	return (x*S_function (x, y, eta, params_ptr));
}

double y_mean(double x, double y, double eta, void *params_ptr)
{
	return (y*S_function (x, y, eta, params_ptr));
}

double eps_3_centershift_real(double x, double y, double eta, void *params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	//note: does not include normalization!!
	double r = sqrt(x*x+y*y);
	double cos_phi = x/r;
	double cos_3_phi = 4.*cos_phi*cos_phi*cos_phi - 3.*cos_phi;

	return (cos_3_phi*S_function (x+params.xavg, y+params.yavg, eta, params_ptr));
}

double eps_3_centershift_imag(double x, double y, double eta, void *params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	//note: does not include normalization!!
	double r = sqrt(x*x+y*y);
	double sin_phi = y/r;
	double sin_3_phi = 3.*sin_phi - 4.*sin_phi*sin_phi*sin_phi;

	return (sin_3_phi*S_function (x+params.xavg, y+params.yavg, eta, params_ptr));
}

//End of file
