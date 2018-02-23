#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_sf_bessel.h>

using namespace std;

#include "main.h"

const int n_tau_pts = 31, n_eta_pts = 31, n_r_pts = 31, n_phi_pts = 31;
const int n_qo_pts = 11, n_qs_pts = 11, n_ql_pts = 11;

vector<double> tau_pts, tau_wts, eta_pts, eta_wts;
vector<double> r_pts, r_wts, phi_pts, phi_wts;
vector<double> S_vector;

vector<double> qo_pts, qs_pts, ql_pts;
double *** correlation_function;
double spectra;

double M = 0.13957, K_Y = 0.0;
double T_0 = 0.12, eta_0 = 0.0, eta_f = 0.6, Delta_eta = 1.2, Rad = 5.0, tau_f = 6.0, Delta_tau = 1.0;
double v_2_bar = 0.0, psi_2_bar = 0.0, eps_2_bar = 0.0;
double v_3_bar = 0.0, psi_3_bar = 0.0, eps_3_bar = 0.0;

double q_max = 0.100;	//GeV
double qo_min = -q_max, qo_max = q_max;
double qs_min = -q_max, qs_max = q_max;
double ql_min = -q_max, ql_max = q_max;

const double tau_min = 0.0, tau_max = 25.0;
const double eta_min = 0.0, eta_max = 4.0;
const double r_min = 0.0, r_max = 50.0;
const double phi_min = 0.0, phi_max = 2*M_PI;

/////////////////////////////////
int main(int argc, char *argv[])
{
	cout << "Starting..." << endl;
	vector<double> KT_pts(21);
	vector<double> KPhi_pts(36);
	linspace(KT_pts, 0.0, 1.0);
	linspace(KPhi_pts, 0.0, 2.0*M_PI);

	#pragma omp parallel for
	for (int iKT = 0; iKT < KT_pts.size(); ++iKT)
	for (int iKphi = 0; iKphi < KPhi_pts.size(); ++iKphi)
	{
		cout << "K_T = " << KT_pts[iKT] << ", K_phi = " << KPhi_pts[iKphi] << endl;
		cout << "\t --> Setting up..." << endl;
		HBT::set_up(M, KT_pts[iKT], KPhi_pts[iKphi], K_Y);

		cout << "\t --> Getting spectra..." << endl;
		HBT::calculate_spectra();

		cout << "\t --> Getting correlation function..." << endl;
		HBT::calculate_correlation_function();

		cout << "\t --> Fitting correlation function..." << endl;
		HBT::fit_correlation_function();

		//output_results();

		cout << "\t --> Cleaning up..." << endl;
		HBT::clean_up();
		cout << endl << endl;
	}

	cout << "Finished all." << endl;
	return 0;
}

// End of file
