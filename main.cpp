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
double ** lambda_CF, ** R2o_CF, ** R2s_CF, ** R2l_CF, ** R2os_CF;
double ** R2o_SV, ** R2s_SV, ** R2l_SV, ** R2os_SV;

void initialize_HBT_data();

/////////////////////////////////
int main(int argc, char *argv[])
{
	cout << "Starting..." << endl;
	vector<double> KT_pts(nKT);
	vector<double> KPhi_pts(nKphi);
	linspace(KT_pts, 0.0, 1.0);
	linspace(KPhi_pts, 0.0, 2.0*M_PI);

	initialize_HBT_data();

	#pragma omp parallel for ordered collapse(2) \
				default(none) shared(M, KT_pts, KPhi_pts, K_Y, correlation_function,\
									R2o_CF, R2s_CF, R2l_CF, R2os_CF)
	for (int iKT = 0; iKT < KT_pts.size(); ++iKT)
	for (int iKphi = 0; iKphi < KPhi_pts.size(); ++iKphi)
	{
		double * results_SV = new double [4];
		double * results_CF = new double [5];
		HBT hbt_corrfunc(M, KT_pts[iKT], KPhi_pts[iKphi], K_Y);

		hbt_corrfunc.calculate_spectra();

		hbt_corrfunc.calculate_SV_radii(results_SV);

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
		cout << KT_pts[iKT] << "   "
				<< KPhi_pts[iKphi] << "   "
				<< K_Y << "   "
				<< qo_pts[iqo] << "   "
				<< qs_pts[iqs] << "   "
				<< ql_pts[iql] << "   "
				<< correlation_function[iKT][iKphi][iqo][iqs][iql] << endl;

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

// End of file
