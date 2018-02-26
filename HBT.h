#ifndef HBT_H
#define HBT_H

using namespace std;

#include "parameters.h"


typedef struct
{
	//double tau, r, phi, eta;
	double t, x, y, z;
	double S_with_weights;
} EmissionFunction;

class HBT
{
	private:
		vector<double> tau_pts, tau_wts, eta_pts, eta_wts;
		vector<double> r_pts, r_wts, phi_pts, phi_wts;
		vector<double> S_vector;

		vector<double> qo_pts, qs_pts, ql_pts;
		double *** correlation_function;
		double spectra;

		double M, K_T, K_Phi, K_Y;
		vector<EmissionFunction> EmissionFunction_vector;

	///////////
	public:
		HBT(double M, double K_T, double K_Phi, double K_Y);

		void calculate_correlation_function(double *** correlation_function);
		void fit_correlation_function(double * results);
		void calculate_spectra();
		void calculate_SV_radii(double * results);
		void clean_up();
		///////////
		double eta_t(double r, double phi);
		double S_function (double tau, double eta, double r, double phi);
		double CF_function(double qt, double qo, double qs, double ql);
};
// End of file

#endif
