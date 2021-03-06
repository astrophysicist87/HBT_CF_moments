#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_errno.h>

#include "gauss/gauss_quadrature.h"
#include "lib.h"

using namespace std;

#define USE_RAPIDITY_SYMMETRY		1

const double hbarC = 0.197327053;		//GeV*fm
const complex<double> complex_i(0.0, 1.0);

const int n_tau_pts = 31, n_eta_pts = 31, n_r_pts = 31, n_phi_pts = 31;
const int n_q_osl_pts = 11;
const int n_qo_pts = n_q_osl_pts, n_qs_pts = n_q_osl_pts, n_ql_pts = n_q_osl_pts;

const double T_0 = 0.12, eta_0 = 0.0, eta_f = 0.6, Delta_eta = 1.2, Rad = 5.0, tau_f = 6.0, Delta_tau = 1.0;
const double v_2_bar = 0.0, psi_2_bar = 0.0, eps_2_bar = 0.0;
const double v_3_bar = 0.1, psi_3_bar = 0.0, eps_3_bar = 0.0;

const double q_max = 0.100;	//GeV
const double qo_min = -q_max, qo_max = q_max;
const double qs_min = -q_max, qs_max = q_max;
const double ql_min = -q_max, ql_max = q_max;

const double tau_min = 0.0, tau_max = 25.0;
const double eta_min = 0.0, eta_max = 4.0;
const double r_min = 0.0, r_max = 50.0;
const double phi_min = 0.0, phi_max = 2.0*M_PI;

// End of file

#endif
