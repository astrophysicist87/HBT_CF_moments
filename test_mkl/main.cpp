#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include "omp.h"
#include "mkl_vml.h"
#include "Stopwatch.h"

using namespace std;

double * v, * result;

/////////////////////////////////
int main(int argc, char *argv[])
{
	int len = 10000000;
	v = new double [len];
	result = new double [len];

	v[0] = 1.0;
	for (int i = 0; i < len-1; i++)
		v[i+1] = 1.0+1.0/(1.0+v[i]);

	Stopwatch sw;
	cout << "Starting main calculation #1..." << endl;
	sw.Start();
	for (int i = 0; i < 500; ++i)
	for (int j = 0; j < len; ++j)
		result[j] = cos(v[j]);
	sw.Stop();
	cout << "Finished main calculation #1 in " << sw.printTime() << " seconds." << endl;

	cout << "Starting main calculation #2..." << endl;
	sw.Start();
	//#pragma omp parallel for
	//for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 500; ++i)
		vdCos (len, v, result);
	sw.Stop();
	cout << "Finished main calculation #2 in " << sw.printTime() << " seconds." << endl;

	delete [] v;
	delete [] result;

	return 0;
}

// End of file
