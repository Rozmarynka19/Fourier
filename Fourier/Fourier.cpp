#define _USE_MATH_DEFINES
#define MaxFactor 13
#include <iostream>
#include <complex>
#include <cmath>
#include <ctime>
using namespace std;
const complex<double> i = 0.0 + 1.0i;

complex<double>* DFT(double* f, int N)
{
	complex<double>* c = new complex<double>[N];
	for (int k = 0; k < N; k++)
	{
		for (int n = 0; n < N; n++)
		{
			c[k] = c[k] + f[n] * exp(-2 * M_PI * i * (double)k * (double)(n /(double)N));
		}
	}
	return c;
}

complex<double>* FFT(double* f, int N)
{
	if (N == 1)
	{
		complex<double>* c = new complex<double>(f[0]);
		return c;
	}
	double* fE = new double[N / 2];
	double* fO = new double[N / 2];
	for (int n = 0; n < N/2; n++)
	{
		fE[n] = f[2 * n];
		fO[n] = f[2 * n + 1];
	}

	complex<double>* cE = FFT(fE, N / 2);
	complex<double>* cO = FFT(fO, N / 2);

	delete[] fE;
	delete[] fO;

	complex<double>* c = new complex<double>[N];

	for (int k = 0; k <= (N / 2 - 1); k++)
	{
		c[k] = cE[k] + exp(-2 * M_PI * i * (double)(k /(double) N)) * cO[k];
	}

	for (int k = 0; k <= (N / 2 - 1); k++)
	{
		c[k+N/2] = cE[k] - exp(-2 * M_PI * i * (double)(k / (double)N)) * cO[k];
	}
	
	delete[] cE;
	delete[] cO;

	return c;
}

double Error(complex<double>* cDFT, complex<double>* cFFT, int N)
{
	double result = 0.0;
	for (int k = 0; k < N; k++)
	{
		result = result + abs(cDFT[k] - cFFT[k]);
	}
	result = result * 1 / (double)N;
	return result;
}

int main()
{
	const bool PrintResult = false;

	for (int i = 1; i <= MaxFactor; i++)
	{
		const int N = 1 << i;
		cout << "N = " << N << endl;

		double* f = new double[N];
		for (int n = 0; n < N; n++)
		{
			f[n] = n /(double)N;
			//f[n] = sin(2 * M_PI * n /N);
		}

		clock_t start = clock();
		complex<double>* cDFT = DFT(f, N);
		clock_t stop = clock();
		double DFTtime = (stop - start) / (double)CLOCKS_PER_SEC;
		cout << "DFT time [s]: " << DFTtime << endl;

		start = clock();
		complex<double>* cFFT = FFT(f, N);
		stop = clock();
		double FFTtime = (stop - start) / (double)CLOCKS_PER_SEC;
		cout << "FFT time [s]: " << FFTtime << endl;

		cout << "Error: " << Error(cDFT, cFFT, N) << endl;

		if (PrintResult)
		{
			cout << "DFT: " << endl;
			for (int n = 0; n < N; n++)
				cout << cDFT[n] << " ";
			cout << "\nFFT: " << endl;
			for (int n = 0; n < N; n++)
				cout << cFFT[n] << " ";
		}

		delete[] f;
		delete[] cDFT;
		delete[] cFFT;
		cout << "\n=========================================\n";
	}
}
//
//int main()
//{
//	//complex<double> z1 = 3.0 + 2.0i;
//	//complex<double> z2 = 5.0 + 3.0i;
//	int N = 8192;
//	double* f = new double[N];
//	cout << "Creating input\n";
//	for (int n = 0; n < N; n++)
//	{
//		f[n] = n / (double)N;
//		//f[n] = sin(2 * M_PI * n /N);
//	}
//	/*for (int n = 0; n < N; n++)
//	{
//		cout << f[n] << " ";
//	}*/
//	cout << "FFT starts...\n";
//	complex<double>* cFFT = DFT(f, N);
//
//	for (int n = 0; n < N; n++)
//	{
//		cout << cFFT[n] << " ";
//	}
//
//	delete[] f;
//	delete[] cFFT;
//
//
//}