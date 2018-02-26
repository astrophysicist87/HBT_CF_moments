#ifndef MAIN_H
#define MAIN_H

#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>

#include "gauss/gauss_quadrature.h"
#include "lib.h"

using namespace std;

#define USE_RAPIDITY_SYMMETRY		1

const double hbarC = 0.197327053;		//GeV*fm

typedef struct
{
	//double tau, r, phi, eta;
	double t, x, y, z;
	double S_with_weights;
} EmissionFunction;

extern const int n_tau_pts, n_eta_pts, n_r_pts, n_phi_pts;
extern const int n_qo_pts, n_qs_pts, n_ql_pts;

extern vector<double> tau_pts, tau_wts, eta_pts, eta_wts;
extern vector<double> r_pts, r_wts, phi_pts, phi_wts;
extern vector<double> qo_pts, qs_pts, ql_pts;
extern vector<double> S_vector;

extern double *** correlation_function;
extern double spectra;

extern double T_0, eta_0, eta_f, Delta_eta, Rad, tau_f, Delta_tau;
extern double v_2_bar, psi_2_bar, eps_2_bar;
extern double v_3_bar, psi_3_bar, eps_3_bar;
extern double qo_min, qo_max, qs_min, qs_max, ql_min, ql_max;
extern const double tau_min, tau_max;
extern const double eta_min, eta_max;
extern const double r_min, r_max;
extern const double phi_min, phi_max;

namespace HBT
{
	double M, K_T, K_Phi, K_Y;
	vector<EmissionFunction> EmissionFunction_vector;

	///////////
	void set_up();
	void calculate_correlation_function();
	void fit_correlation_function();
	void output_results();
	void clean_up();
	///////////
	double eta_t(double r, double phi);
	double S_function (double tau, double eta, double r, double phi);
	///////////

	void set_up(double M_in, double K_T_in, double K_Phi_in, double K_Y_in)
	{
		M = M_in;
		K_T = K_T_in;
		K_Phi = K_Phi_in;
		K_Y = K_Y_in;
		S_vector.clear();
		EmissionFunction_vector.clear();

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
		{
			double tau = tau_pts[itau];
			double eta = eta_pts[ieta];
			double r = r_pts[ir];
			double phi = phi_pts[iphi];

			double t = tau * cosh(eta);
			double x = r * cos(phi);
			double y = r * sin(phi);
			double z = tau * sinh(eta);

			EmissionFunction local_EF;

			//EmissionFunction_vector[idx].tau = tau;
			//EmissionFunction_vector[idx].r = r;
			//EmissionFunction_vector[idx].phi = phi;
			//EmissionFunction_vector[idx].eta = eta;
			local_EF.t = t;
			local_EF.x = x;
			local_EF.y = y;
			local_EF.z = z;

			S_vector[idx] = S_function(tau, eta, r, phi);
			//note inclusion of tau and r factors in weights
			local_EF.S_with_weights
				= tau_wts[itau]*tau_pts[itau]*eta_wts[ieta]*r_wts[ir]
					*phi_wts[iphi]*r_pts[ir]*S_vector[idx];

			EmissionFunction_vector.push_back(local_EF);

			++idx;
		}

		//set up CF
		correlation_function = new double ** [n_qo_pts];
		for (int iqo = 0; iqo < n_qo_pts; ++iqo)
		{
			correlation_function[iqo] = new double * [n_qs_pts];
			for (int iqs = 0; iqs < n_qs_pts; ++iqs)
			{
				correlation_function[iqo][iqs] = new double [n_ql_pts];
				for (int iql = 0; iql < n_ql_pts; ++iql)
					correlation_function[iqo][iqs][iql] = 0.0;
			}
		
		}
	
		return;
	}

	double eta_t(double r, double phi)
	{
		return (
			eta_f * (r/Rad)
					* ( 1.0
						+ 2.0 * v_2_bar * ( cos( 2.0 * ( phi - psi_2_bar ) ) )
						+ 2.0 * v_3_bar * ( cos( 3.0 * ( phi - psi_3_bar ) ) ) )
				);
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

	void calculate_spectra()
	{
		double sym_factor = ( USE_RAPIDITY_SYMMETRY ) ? 2.0 : 1.0;
		spectra = 0.0;
		int idx = 0;
		for (int itau = 0; itau < n_tau_pts; ++itau)
		for (int ieta = 0; ieta < n_eta_pts; ++ieta)
		for (int ir = 0; ir < n_r_pts; ++ir)
		for (int iphi = 0; iphi < n_phi_pts; ++iphi)
			spectra += sym_factor * eta_wts[ieta]
						* tau_wts[itau] * tau_pts[itau]
						* r_wts[ir] * r_pts[ir] * phi_wts[iphi]
						* S_vector[idx++];
	}

	double CF_function(double qt, double qo, double qs, double ql)
	{
		if (not USE_RAPIDITY_SYMMETRY)
		{
			cerr << "Need to get this function working for this option!!!  Exiting..." << endl;
			exit(8);
		}

		//cout << "Doing (qo, qs, ql) = (" << qo << ", " << qs << ", " << ql << ")..." << endl;
		complex<double> result = 0.0;
		//vector<complex<double> > phase_vector(S_vector.size());

		const int Ncells = (int)S_vector.size();
		double cos_K_Phi = cos(K_Phi), sin_K_Phi = sin(K_Phi);

		//#pragma omp simd
		for (int idx = 0; idx < Ncells; ++idx)
		{
			EmissionFunction * local_EF = &(EmissionFunction_vector[idx]);
			double t = local_EF->t;
			double x = local_EF->x;
			double y = local_EF->y;
			double z = local_EF->z;
			double xo = x*cos_K_Phi+y*sin_K_Phi;
			double xs = y*cos_K_Phi-x*sin_K_Phi;
			double xl = z;

			complex<double> phase
				= ( USE_RAPIDITY_SYMMETRY ) ? 
					2.0*exp( i*( qt*t - ( qo*xo+qs*xs ) )/hbarC ) * cos(ql*xl/hbarC):
					exp( i*( qt*t - ( qo*xo+qs*xs+ql*xl ) )/hbarC );

			result += phase*(local_EF->S_with_weights);
		}
	
		return ( 1.0 + norm( result / spectra ) );
	}

	void calculate_correlation_function()
	{
		#pragma omp parallel for collapse(3)
		for (int iqo = 0; iqo < n_qo_pts; ++iqo)
		for (int iqs = 0; iqs < n_qs_pts; ++iqs)
		for (int iql = 0; iql < n_ql_pts; ++iql)
		{
			double qo = qo_pts[iqo], qs = qs_pts[iqs], ql = ql_pts[iql];
			double q_dot_K = qo*K_T;
			double xi2 = M*M+K_T*K_T+0.25*(qo*qo+qs*qs+ql*ql);
			double E1 = sqrt(xi2 + q_dot_K);
			double E2 = sqrt(xi2 - q_dot_K);
			double qt = E1 - E2;
			correlation_function[iqo][iqs][iql] = CF_function(qt, qo, qs, ql);
		}

		return;
	}

	void fit_correlation_function()
	{
		const size_t data_length = n_qo_pts*n_qs_pts*n_ql_pts;  // # of points

		double lambda, R_o, R_s, R_l, R_os;
		int dim = 5;
		int s_gsl;

		double * V = new double [dim];
		double * qweight = new double [dim];
		double ** T = new double * [dim];
		for(int i = 0; i < dim; i++)
		{
		    V[i] = 0.0;
		    T[i] = new double [dim];
		    for(int j = 0; j < dim; j++)
		        T[i][j] = 0.0;
		}

		gsl_matrix * T_gsl = gsl_matrix_alloc (dim, dim);
		gsl_matrix * T_inverse_gsl = gsl_matrix_alloc (dim, dim);
		gsl_permutation * perm = gsl_permutation_alloc (dim);

		//double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
		double CF_err = 1.e-3;
		for (int i = 0; i < n_qo_pts; i++)
		for (int j = 0; j < n_qs_pts; j++)
		for (int k = 0; k < n_ql_pts; k++)
		{
		    //double qo = q1pts[i] * ckp + q2pts[j] * skp;
		    //double qs = -q1pts[i] * skp + q2pts[j] * ckp;
		    //double ql = q3pts[k];
		    double qo = qo_pts[i];
		    double qs = qs_pts[j];
		    double ql = ql_pts[k];
		    double correl_local = correlation_function[i][j][k]-1;
		    if(correl_local < 1e-15) continue;
		    double sigma_k_prime = CF_err/correl_local;
		        
		    double inv_sigma_k_prime_sq = 1./(sigma_k_prime*sigma_k_prime);
		    double log_correl_over_sigma_sq = log(correl_local)*inv_sigma_k_prime_sq;

		    qweight[0] = - 1.0;
		    qweight[1] = qo*qo;
		    qweight[2] = qs*qs;
		    qweight[3] = ql*ql;
		    qweight[4] = qo*qs;

		    for(int ij = 0; ij < dim; ij++)
		    {
		        V[ij] += qweight[ij]*log_correl_over_sigma_sq;
		        T[0][ij] += qweight[ij]*inv_sigma_k_prime_sq;
		    }

		    for(int ij = 1; ij < dim; ij++)
		        T[ij][0] = T[0][ij];
		        

		    for(int ij = 1; ij < dim; ij++)
		    {
		        for(int lm = 1; lm < dim; lm++)
		            T[ij][lm] += -qweight[ij]*qweight[lm]*inv_sigma_k_prime_sq;
		    }
		}
		for(int i = 0; i < dim; i++)
		    for(int j = 0; j < dim; j++)
		        gsl_matrix_set(T_gsl, i, j, T[i][j]);

		// Make LU decomposition of matrix T_gsl
		gsl_linalg_LU_decomp (T_gsl, perm, &s_gsl);
		// Invert the matrix m
		gsl_linalg_LU_invert (T_gsl, perm, T_inverse_gsl);

		double **T_inverse = new double* [dim];
		for(int i = 0; i < dim; i++)
		{
		    T_inverse[i] = new double [dim];
		    for(int j = 0; j < dim; j++)
		        T_inverse[i][j] = gsl_matrix_get(T_inverse_gsl, i, j);
		}
		double *results = new double [dim];
		for(int i = 0; i < dim; i++)
		{
		    results[i] = 0.0;
		    for(int j = 0; j < dim; j++)
		        results[i] += T_inverse[i][j]*V[j];
		}

	/*
		lambda_Correl[ipt][ipphi] = exp(results[0]);
		lambda_Correl_err[ipt][ipphi] = 0.0;
		R2_out_GF[ipt][ipphi] = results[1]*hbarC*hbarC;
		R2_side_GF[ipt][ipphi] = results[2]*hbarC*hbarC;
		R2_long_GF[ipt][ipphi] = results[3]*hbarC*hbarC;
		R2_outside_GF[ipt][ipphi] = results[4]*hbarC*hbarC;
		R2_out_err[ipt][ipphi] = 0.0;
		R2_side_err[ipt][ipphi] = 0.0;
		R2_long_err[ipt][ipphi] = 0.0;
		R2_outside_err[ipt][ipphi] = 0.0;
	*/
	cout << "lambda = " << exp(results[0]) << endl;
	cout << "R2o = " << results[1]*hbarC*hbarC << endl;
	cout << "R2s = " << results[2]*hbarC*hbarC << endl;
	cout << "R2l = " << results[3]*hbarC*hbarC << endl;
	cout << "R2os = " << results[4]*hbarC*hbarC << endl;
	cout << "RESULTS: "
			<< K_T << "   " << K_Phi << "   " << K_Y << "   "
			<< exp(results[0]) << "   "
			<< results[1]*hbarC*hbarC << "   "
			<< results[2]*hbarC*hbarC << "   "
			<< results[3]*hbarC*hbarC << "   "
			<< results[4]*hbarC*hbarC << endl;

		double chi_sq = 0.0;
		for (int i = 0; i < n_qo_pts; i++)
		for (int j = 0; j < n_qs_pts; j++)
		for (int k = 0; k < n_ql_pts; k++)
		{
		    double qo = qo_pts[i];
		    double qs = qs_pts[j];
		    double ql = ql_pts[k];
		    double correl_local = correlation_function[i][j][k]-1;
		    if(correl_local < 1e-15) continue;
		    double sigma_k_prime = CF_err/correl_local;

		    chi_sq += pow((log(correl_local) - results[0] 
		                   + results[1]*qo*qo 
		                   + results[2]*qs*qs
		                   + results[3]*ql*ql
		                   + results[4]*qo*qs), 2)
		              /sigma_k_prime/sigma_k_prime;
		}
		//cout << "chi_sq/d.o.f = " << chi_sq/(qnpts - dim) << endl;
		//chi_sq_per_dof = chi_sq/(qnpts - dim);

		// clean up
		gsl_matrix_free (T_gsl);
		gsl_matrix_free (T_inverse_gsl);
		gsl_permutation_free (perm);

		delete [] qweight;
		delete [] V;
		for(int i = 0; i < dim; i++)
		{
		    delete [] T[i];
		    delete [] T_inverse[i];
		}
		delete [] T;
		delete [] T_inverse;
		delete [] results;
	}

	void clean_up()
	{
		//clean up CF
		for (int iqo = 0; iqo < n_qo_pts; ++iqo)
		{
			for (int iqs = 0; iqs < n_qs_pts; ++iqs)
				delete [] correlation_function[iqo][iqs];
			delete [] correlation_function[iqo];
		}
		delete [] correlation_function;

	}

}
// End of file

#endif
