#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

using namespace std;

#define USE_RAPIDITY_SYMMETRY		1

extern const int n_tau_pts, n_eta_pts, n_r_pts, n_phi_pts;
extern const int n_qo_pts, n_qs_pts, n_ql_pts;

extern vector<double> tau_pts, tau_wts, eta_pts, eta_wts;
extern vector<double> r_pts, r_wts, phi_pts, phi_wts;
extern vector<double> qo_pts, qs_pts, ql_pts;

extern double M, K_Y, K_T, K_Phi;
extern double T_0, eta_0, Delta_eta, Rad, tau_f, Delta_tau;
extern double v_2_bar, psi_2_bar, eps_2_bar;
extern double v_3_bar, psi_3_bar, eps_3_bar;
extern double qo_min, qo_max, qs_min, qs_max, ql_min, ql_max;

void set_up()
{
	S_vector = vector<double>(n_tau_pts*n_eta_pts*n_r_pts*n_phi_pts);

	tau_pts = vector<double>(n_tau_pts);
	tau_wts = vector<double>(n_tau_pts);
	gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, tau_min, tau_max, tau_pts, tau_wts);

	eta_pts = vector<double>(n_eta_pts);
	eta_wts = vector<double>(n_eta_pts);
	gauss_quadrature(n_eta_pts, 1, 0.0, 0.0, eta_min, eta_max, eta_pts, eta_wts);

	r_pts = vector<double>(n_r_pts);
	r_wts = vector<double>(n_r_pts);
	gauss_quadrature(n_r_pts, 1, 0.0, 0.0, r_min, r_max, r_pts, r_wts);

	phi_pts = vector<double>(n_phi_pts);
	phi_wts = vector<double>(n_phi_pts);
	gauss_quadrature(n_phi_pts, 1, 0.0, 0.0, phi_min, phi_max, phi_pts, phi_wts);

	qo_pts = vector<double>(n_qo_pts);
	qs_pts = vector<double>(n_qs_pts);
	ql_pts = vector<double>(n_ql_pts);
	linspace(qo_pts, qo_min, qo_max);
	linspace(qs_pts, qs_min, qs_max);
	linspace(ql_pts, ql_min, ql_max);

	int idx = 0;
	for (int itau = 0; itau < n_tau_pts; ++itau)
	for (int ieta = 0; ieta < n_eta_pts; ++ieta)
	for (int ir = 0; ir < n_r_pts; ++ir)
	for (int iphi = 0; iphi < n_phi_pts; ++iphi)
		S_vector[idx++] = S_function(tau_pts[itau], eta_pts[ieta], r_pts[ir], phi_pts[iphi]);

	return;
}

double S_function (double tau, double eta, double r, double phi)
{
	double M_T = sqrt(M*M+K_T*K_T);
	double eta_t_loc = eta_t(r, phi);

	//Expression for function is long and complicated,
	//so break it up into constituent terms
	double term0, term1, term2, term3, term4;
	term0 = (tau-tau_f)*(tau-tau_f) / (2.*Delta_tau*Delta_tau);
	term1 = (eta-eta_0)*(eta-eta_0) / (2.*Delta_eta*Delta_eta);
	term2 = (M_T/T_0) * cosh(eta - K_Y) * cosh(eta_t_loc);
	term3 = (K_T/T_0) * (cos(phi - K_Phi)) * sinh(eta_t_loc);
	term4 = (r*r)/(2.*Rad*Rad)
				* ( 1.0
					+ 2.0 * eps_2_bar * ( cos( 2.0 * ( phi - psi_2_bar ) ) )
					+ 2.0 * eps_3_bar * ( cos( 3.0 * ( phi - psi_3_bar ) ) )
					);

	return (exp(-term0 - term1 - term2 + term3 - term4));
}


double eta_t(double r, double phi)
{
	return (
		eta_f * (r/Rad)
				* ( 1.0
					+ 2.0 * v_2_bar * ( cos( 2.0 * ( phi - psi_2_bar ) ) )
					+ 2.0 * v_3_bar * ( cos( 3.0 * ( phi - psi_3_bar ) ) )
			);
}

double CF_function(double qt, double qo, double qs, double ql)
{
	if (not USE_RAPIDITY_SYMMETRY)
	{
		cerr << "Need to get this function working for this option!!!  Exiting..." << endl;
		exit(8);
	}

	int idx = 0;
	complex<double> result = 0.0;
	for (int itau = 0; itau < n_tau_pts; ++itau)
	for (int ieta = 0; ieta < n_eta_pts; ++ieta)
	for (int ir = 0; ir < n_r_pts; ++ir)
	for (int iphi = 0; iphi < n_phi_pts; ++iphi)
	{
		double tau = tau_pts[itau];
		double eta = eta_pts[ieta];
		double r = r_pts[ir];
		double phi = phi_pts[iphi];

		double t = tau * cosh(eta);
		double xo = r * cos(phi - K_Phi);
		double xs = r * sin(phi - K_Phi);
		double xl = tau * sinh(eta);

		complex<double> phase
			= ( USE_RAPIDITY_SYMMETRY ) ? 
				2.0*exp( i*( qt*t - ( qo*xo+qs*xs ) ) ) * cos(ql*xl):
				exp( i*( qt*t - ( qo*xo+qs*xs+ql*xl ) ) );

		result += tau_wts[itau] * eta_wts[ieta]
					* r_wts[ir] * phi_wts[iphi]
					* phase * S_vector[idx++];
	}
	
	return ( 1.0 + norm( result ) );
}

void calculate_correlation_function()
{
	for (int iqo = 0; iqo < n_qo_pts; ++iqo)
	for (int iqs = 0; iqs < n_qs_pts; ++iqs)
	for (int iql = 0; iql < n_ql_pts; ++iql)
	{
		double E1 = ...;
		double E2 = ...;
		double qt = E2 - E1;
		correlation_function[iqo][iqs][iql] = CF_function(qt, qo_pts[iqo], qs_pts[iqs], ql_pts[iql]);
	}

	return;
}

void fit_correlation_function()
{
	...;
}

// End of file

#endif
