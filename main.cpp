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

double M = 0.13957, K_Y = 0.0;

/////////////////////////////////
int main(int argc, char *argv[])
{
	cout << "Starting..." << endl;
	vector<double> KT_pts(21);
	vector<double> KPhi_pts(36);
	linspace(KT_pts, 0.0, 1.0);
	linspace(KPhi_pts, 0.0, 2.0*M_PI);

	#pragma omp parallel for ordered collapse(2) \
				default(none) shared(M, KT_pts, KPhi_pts, K_Y)
	for (int iKT = 0; iKT < KT_pts.size(); ++iKT)
	for (int iKphi = 0; iKphi < KPhi_pts.size(); ++iKphi)
	{

		HBT hbt_corrfunc(M, KT_pts[iKT], KPhi_pts[iKphi], K_Y);

		hbt_corrfunc.calculate_spectra();

		hbt_corrfunc.calculate_correlation_function();

		hbt_corrfunc.fit_correlation_function();

		//output_results();

		hbt_corrfunc.clean_up();

	}

	cout << "Finished all." << endl;
	return 0;
}

// End of file
