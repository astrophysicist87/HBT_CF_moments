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

#include "parameters.h"
#include "lib.h"
#include "HBT.h"

const int nKT = 3, nKphi = 4;
double M = 0.13957, K_Y = 0.0;
double ***** correlation_function;
vector<double> KT_pts;
vector<double> KPhi_pts, KPhi_wts;
double ** lambda_CF, ** R2o_CF, ** R2s_CF, ** R2l_CF, ** R2os_CF;
double ** R2o_CF_C, ** R2s_CF_C, ** R2l_CF_C, ** R2os_CF_C;
double ** R2o_CF_S, ** R2s_CF_S, ** R2l_CF_S, ** R2os_CF_S;
double ** R2o_SV, ** R2s_SV, ** R2l_SV, ** R2os_SV;
double ** R2o_SV_C, ** R2s_SV_C, ** R2l_SV_C, ** R2os_SV_C;
double ** R2o_SV_S, ** R2s_SV_S, ** R2l_SV_S, ** R2os_SV_S;

void initialize_HBT_data();
void R2ij_Fourier_transform(int n_order, int SV_or_CF_mode);

/////////////////////////////////
int main(int argc, char *argv[])
{
	cout << "Starting..." << endl;
	KT_pts = vector<double>(nKT);
	KPhi_pts = vector<double>(nKphi);
	KPhi_wts = vector<double>(nKphi);
	linspace(KT_pts, 0.0, 1.0);
	//linspace(KPhi_pts, 0.0, 2.0*M_PI);
	gauss_quadrature(nKphi, 1, 0.0, 0.0, 0.0, 2.0*M_PI, KPhi_pts, KPhi_wts);

	initialize_HBT_data();

	#pragma omp parallel for ordered collapse(2) \
				default(none) shared(M, KT_pts, KPhi_pts, K_Y, correlation_function,\
									lambda_CF, R2o_CF, R2s_CF, R2l_CF, R2os_CF,\
									R2o_SV, R2s_SV, R2l_SV, R2os_SV)
	for (int iKT = 0; iKT < KT_pts.size(); ++iKT)
	for (int iKphi = 0; iKphi < KPhi_pts.size(); ++iKphi)
	{
		double * results_SV = new double [4];
		double * results_CF = new double [5];
		HBT hbt_corrfunc(M, KT_pts[iKT], KPhi_pts[iKphi], K_Y);

		hbt_corrfunc.calculate_spectra();

		hbt_corrfunc.calculate_SV_radii(results_SV);
		R2o_SV[iKT][iKphi] = results_SV[0];
		R2s_SV[iKT][iKphi] = results_SV[1];
		R2l_SV[iKT][iKphi] = results_SV[2];
		R2os_SV[iKT][iKphi] = results_SV[3];

		hbt_corrfunc.calculate_correlation_function(correlation_function[iKT][iKphi]);

		hbt_corrfunc.fit_correlation_function(results_CF);

		lambda_CF[iKT][iKphi] = results_CF[0];
		R2o_CF[iKT][iKphi] = results_CF[1];
		R2s_CF[iKT][iKphi] = results_CF[2];
		R2l_CF[iKT][iKphi] = results_CF[3];
		R2os_CF[iKT][iKphi] = results_CF[4];

		hbt_corrfunc.clean_up();

		delete [] results_SV;
		delete [] results_CF;
	}

	//
	vector<double> qo_pts(n_qo_pts);
	vector<double> qs_pts(n_qs_pts);
	vector<double> ql_pts(n_ql_pts);
	linspace(qo_pts, qo_min, qo_max);
	linspace(qs_pts, qs_min, qs_max);
	linspace(ql_pts, ql_min, ql_max);

	//
	for (int iKT = 0; iKT < KT_pts.size(); ++iKT)
	for (int iKphi = 0; iKphi < KPhi_pts.size(); ++iKphi)
	for (int iqo = 0; iqo < n_qo_pts; ++iqo)
	for (int iqs = 0; iqs < n_qs_pts; ++iqs)
	for (int iql = 0; iql < n_ql_pts; ++iql)
		cout << "CF-1: " << KT_pts[iKT] << "   "
				<< KPhi_pts[iKphi] << "   "
				<< K_Y << "   "
				<< qo_pts[iqo] << "   "
				<< qs_pts[iqs] << "   "
				<< ql_pts[iql] << "   "
				<< correlation_function[iKT][iKphi][iqo][iqs][iql] << endl;

	cout << endl;
	//
	for (int iKT = 0; iKT < KT_pts.size(); ++iKT)
	for (int iKphi = 0; iKphi < KPhi_pts.size(); ++iKphi)
		cout << "R2ij(SV): " << KT_pts[iKT] << "   "
				<< KPhi_pts[iKphi] << "   "
				<< K_Y << "   "
				<< R2o_SV[iKT][iKphi] << "   "
				<< R2s_SV[iKT][iKphi] << "   "
				<< R2l_SV[iKT][iKphi] << "   "
				<< R2os_SV[iKT][iKphi] << endl;

	cout << endl;
	//
	for (int iKT = 0; iKT < KT_pts.size(); ++iKT)
	for (int iKphi = 0; iKphi < KPhi_pts.size(); ++iKphi)
		cout << "R2ij(CF): " << KT_pts[iKT] << "   "
				<< KPhi_pts[iKphi] << "   "
				<< K_Y << "   "
				<< R2o_CF[iKT][iKphi] << "   "
				<< R2s_CF[iKT][iKphi] << "   "
				<< R2l_CF[iKT][iKphi] << "   "
				<< R2os_CF[iKT][iKphi] << endl;


	cout << "Finished all." << endl;
	return 0;
}

void initialize_HBT_data()
{
	correlation_function = new double **** [nKT];
	lambda_CF = new double * [nKT];
	R2o_CF = new double * [nKT];
	R2s_CF = new double * [nKT];
	R2l_CF = new double * [nKT];
	R2os_CF = new double * [nKT];
	R2o_SV = new double * [nKT];
	R2s_SV = new double * [nKT];
	R2l_SV = new double * [nKT];
	R2os_SV = new double * [nKT];
	for (int iKT = 0; iKT < nKT; ++iKT)
	{
		correlation_function[iKT] = new double *** [nKphi];
		lambda_CF[iKT] = new double [nKphi];
		R2o_CF[iKT] = new double [nKphi];
		R2s_CF[iKT] = new double [nKphi];
		R2l_CF[iKT] = new double [nKphi];
		R2os_CF[iKT] = new double [nKphi];
		R2o_SV[iKT] = new double [nKphi];
		R2s_SV[iKT] = new double [nKphi];
		R2l_SV[iKT] = new double [nKphi];
		R2os_SV[iKT] = new double [nKphi];
		for (int iKphi = 0; iKphi < nKphi; ++iKphi)
		{
			correlation_function[iKT][iKphi] = new double ** [n_qo_pts];
			lambda_CF[iKT][iKphi] = 0.0;
			R2o_CF[iKT][iKphi] = 0.0;
			R2s_CF[iKT][iKphi] = 0.0;
			R2l_CF[iKT][iKphi] = 0.0;
			R2os_CF[iKT][iKphi] = 0.0;
			R2o_SV[iKT][iKphi] = 0.0;
			R2s_SV[iKT][iKphi] = 0.0;
			R2l_SV[iKT][iKphi] = 0.0;
			R2os_SV[iKT][iKphi] = 0.0;
			for (int iqo = 0; iqo < n_qo_pts; ++iqo)
			{
				correlation_function[iKT][iKphi][iqo] = new double * [n_qs_pts];
				for (int iqs = 0; iqs < n_qs_pts; ++iqs)
				{
					correlation_function[iKT][iKphi][iqo][iqs] = new double [n_ql_pts];
					for (int iql = 0; iql < n_ql_pts; ++iql)
						correlation_function[iKT][iKphi][iqo][iqs][iql] = 0.0;
				}
			}
		}
	}
	return;
}

void R2ij_Fourier_transform(int n_order, int SV_or_CF_mode)
{
	double ** R2s_mode, ** R2o_mode, ** R2l_mode, ** R2os_mode;
	double ** R2s_mode_C, ** R2o_mode_C, ** R2l_mode_C, ** R2os_mode_C;
	double ** R2s_mode_S, ** R2o_mode_S, ** R2l_mode_S, ** R2os_mode_S;
	switch (SV_or_CF_mode)
	{
		case 0:
			R2s_mode = R2s_SV;
			R2o_mode = R2o_SV;
			R2l_mode = R2l_SV;
			R2os_mode = R2os_SV;
			R2s_mode_C = R2s_SV_C;
			R2o_mode_C = R2o_SV_C;
			R2l_mode_C = R2l_SV_C;
			R2os_mode_C = R2os_SV_C;
			R2s_mode_S = R2s_SV_S;
			R2o_mode_S = R2o_SV_S;
			R2l_mode_S = R2l_SV_S;
			R2os_mode_S = R2os_SV_S;
			break;
		case 1:
			R2s_mode = R2s_CF;
			R2o_mode = R2o_CF;
			R2l_mode = R2l_CF;
			R2os_mode = R2os_CF;
			R2s_mode_C = R2s_CF_C;
			R2o_mode_C = R2o_CF_C;
			R2l_mode_C = R2l_CF_C;
			R2os_mode_C = R2os_CF_C;
			R2s_mode_S = R2s_CF_S;
			R2o_mode_S = R2o_CF_S;
			R2l_mode_S = R2l_CF_S;
			R2os_mode_S = R2os_CF_S;
			break;
		default:
			cerr << "Not a supported value!" << endl;
			exit (1);
			break;
	}

	for(int Morder = 0; Morder < n_order; ++Morder)
	{
		//set plane_psi...
		double plane_psi = 0.0;	//for now

		double cos_mKPhi_pts[nKphi], sin_mKPhi_pts[nKphi];

		for(int iKphi = 0; iKphi < nKphi; ++iKphi)
		{
			cos_mKPhi_pts[iKphi] = cos(Morder*(KPhi_pts[iKphi] - plane_psi));
			sin_mKPhi_pts[iKphi] = sin(Morder*(KPhi_pts[iKphi] - plane_psi));
		}

		for(int iKT = 0; iKT < nKT; ++iKT)
		{
			double temp_sum_side_cos = 0.0, temp_sum_side_sin = 0.0;
			double temp_sum_out_cos = 0.0, temp_sum_out_sin = 0.0;
			double temp_sum_outside_cos = 0.0, temp_sum_outside_sin = 0.0;
			double temp_sum_long_cos = 0.0, temp_sum_long_sin = 0.0;

			for(int iKphi = 0; iKphi < nKphi; ++iKphi)
			{
				double local_R2s = R2s_mode[iKT][iKphi];
				double local_R2o = R2o_mode[iKT][iKphi];
				double local_R2l = R2l_mode[iKT][iKphi];
				double local_R2os = R2os_mode[iKT][iKphi];
			
				temp_sum_side_cos += local_R2s*cos_mKPhi_pts[iKphi]*KPhi_wts[iKphi];
				temp_sum_side_sin += local_R2s*sin_mKPhi_pts[iKphi]*KPhi_wts[iKphi];
				temp_sum_out_cos += local_R2o*cos_mKPhi_pts[iKphi]*KPhi_wts[iKphi];
				temp_sum_out_sin += local_R2o*sin_mKPhi_pts[iKphi]*KPhi_wts[iKphi];
				temp_sum_outside_cos += local_R2os*cos_mKPhi_pts[iKphi]*KPhi_wts[iKphi];
				temp_sum_outside_sin += local_R2os*sin_mKPhi_pts[iKphi]*KPhi_wts[iKphi];
				temp_sum_long_cos += local_R2l*cos_mKPhi_pts[iKphi]*KPhi_wts[iKphi];
				temp_sum_long_sin += local_R2l*sin_mKPhi_pts[iKphi]*KPhi_wts[iKphi];
			}

			R2s_mode_C[iKT][Morder] = temp_sum_side_cos/(2.*M_PI);
			R2s_mode_S[iKT][Morder] = temp_sum_side_sin/(2.*M_PI);
			R2o_mode_C[iKT][Morder] = temp_sum_out_cos/(2.*M_PI);
			R2o_mode_S[iKT][Morder] = temp_sum_out_sin/(2.*M_PI);
			R2os_mode_C[iKT][Morder] = temp_sum_outside_cos/(2.*M_PI);
			R2os_mode_S[iKT][Morder] = temp_sum_outside_sin/(2.*M_PI);
			R2l_mode_C[iKT][Morder] = temp_sum_long_cos/(2.*M_PI);
			R2l_mode_S[iKT][Morder] = temp_sum_long_sin/(2.*M_PI);
		}
	}

	return;

}

// End of file
