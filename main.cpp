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

using namespace std;

#include "parameters.h"
#include "lib.h"
#include "HBT_at_K.h"
#include "main.h"

/////////////////////////////////
int main(int argc, char *argv[])
{
	gsl_set_error_handler_off();

	cout << "Starting..." << endl;
	KT_pts = vector<double>(nKT);
	KPhi_pts = vector<double>(nKphi);
	KPhi_wts = vector<double>(nKphi);
	linspace(KT_pts, 0.0, 1.0);
	//linspace(KPhi_pts, 0.0, 2.0*M_PI);
	gauss_quadrature(nKphi, 1, 0.0, 0.0, 0.0, 2.0*M_PI, KPhi_pts, KPhi_wts);

	initialize_HBT_data();

	//begin parallel
	#pragma omp parallel for ordered collapse(2) \
				default(none) shared(M, KT_pts, KPhi_pts, K_Y, correlation_function,\
									lambda_CF, R2o_CF, R2s_CF, R2l_CF, R2os_CF,\
									R2o_SV, R2s_SV, R2l_SV, R2os_SV, spectra)
	for (int iKT = 0; iKT < KT_pts.size(); ++iKT)
	for (int iKphi = 0; iKphi < KPhi_pts.size(); ++iKphi)
	{
		double * results_SV = new double [4];
		double * results_CF = new double [5];
		HBT_at_K hbt_at_k_corrfunc(M, KT_pts[iKT], KPhi_pts[iKphi], K_Y);

		hbt_at_k_corrfunc.calculate_spectra(spectra[iKT][iKphi]);

		if (COMPUTE_HBT_SV)
		{
			hbt_at_k_corrfunc.calculate_SV_radii(results_SV);
			R2o_SV[iKT][iKphi] = results_SV[0];
			R2s_SV[iKT][iKphi] = results_SV[1];
			R2l_SV[iKT][iKphi] = results_SV[2];
			R2os_SV[iKT][iKphi] = results_SV[3];
		}

		if (COMPUTE_HBT_CF)
		{
			hbt_at_k_corrfunc.calculate_correlation_function(correlation_function[iKT][iKphi]);

			hbt_at_k_corrfunc.fit_correlation_function(results_CF);

			lambda_CF[iKT][iKphi] = results_CF[0];
			R2o_CF[iKT][iKphi] = results_CF[1];
			R2s_CF[iKT][iKphi] = results_CF[2];
			R2l_CF[iKT][iKphi] = results_CF[3];
			R2os_CF[iKT][iKphi] = results_CF[4];
		}

		hbt_at_k_corrfunc.clean_up();

		delete [] results_SV;
		delete [] results_CF;
	}
	//end parallel

	//set flow plane angles
	get_flow();

	//compute Fourier transform of differential radii
	R2ij_Fourier_transform(n_max, 0);
	R2ij_Fourier_transform(n_max, 1);

	//get correlation function Fourier moments
	get_CF_moments();

	//begin parallel
	#pragma omp parallel for ordered collapse(2) \
				default(none) shared(KT_pts, correlation_function_moments,\
									R2o_CFM_C, R2s_CFM_C, R2l_CFM_C, R2os_CFM_C,\
									R2o_CFM_S, R2s_CFM_S, R2l_CFM_S, R2os_CFM_S)
	for (int iKT = 0; iKT < KT_pts.size(); ++iKT)
	for (int order = 0; order <= n_max; ++order)
	{
		double * results_C = new double [5];
		double * results_S = new double [5];

		fit_CF_moments(correlation_function_moments[iKT][order][0], results_C);

		//lambda_CF[iKT][order] = results_C[0];
		R2o_CFM_C[iKT][order] = results_C[1];
		R2s_CFM_C[iKT][order] = results_C[2];
		R2l_CFM_C[iKT][order] = results_C[3];
		R2os_CFM_C[iKT][order] = results_C[4];

		fit_CF_moments(correlation_function_moments[iKT][order][1], results_S);

		//lambda_CF[iKT][order] = results_S[0];
		R2o_CFM_S[iKT][order] = results_S[1];
		R2s_CFM_S[iKT][order] = results_S[2];
		R2l_CFM_S[iKT][order] = results_S[3];
		R2os_CFM_S[iKT][order] = results_S[4];

		delete [] results_C;
		delete [] results_S;
	}
	//end parallel

	print_some_stuff();

	cout << "Finished all." << endl;
	return 0;
}

// End of file
