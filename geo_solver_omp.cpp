#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <stdio.h>
#include <omp.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define N 902
#define R 6378000.0
#define GM 398600.5
#define D 300000//o kolko je vnoreny bod

int main() {
	double r_ij, rx, ry, rz,sum;
	int i, j;

	double** A;
	double* B = new double[N]();
	double* L = new double[N]();
	double* H = new double[N]();
	double* g = new double[N]();
	double* u2n2 = new double[N]();
	double* G_temp = new double[N];//riadok matice na nasobenie
	double* Apj = new double[N]();
	double* As = new double[N]();

	double* X = new double[N];
	double* Y = new double[N];
	double* Z = new double[N];
	double* x = new double[N];
	double* y = new double[N];
	double* z = new double[N];
	double* nx = new double[N];
	double* ny = new double[N];
	double* nz = new double[N];
	double* u = new double[N];

	double* xj;
	double* r0;
	double* r_prev;
	double* pj;
	double* x_prev;
	double* rozdiel;
	double* s;
	double* rj;
	double* temp_result;
	double* temp_result2;

	int iter = 0, maxIterations = 100;
	double alfa = 0, omega = 0, beta = 0, epsilon = 0.00001;




	FILE* file = fopen("BL-8102.dat", "r");

	if (!file) {
		perror("Error opening file");
	}

	for (i = 0; i < N; i++)
	{
		int count = fscanf(file, "%lf %lf %lf %lf %lf", &B[i], &L[i], &H[i], &g[i], &u2n2[i]);
		//g[i] = GM / (R * R);
		//H[i] = 0;
		g[i] = g[i] * 0.00001;


	}
	fclose(file);

	A = new double* [N];
	for ( i = 0; i < N; i++)
	{
		A[i] = new double[N]();
	}

	for ( i = 0; i < N; i++)
	{
		double Brad = B[i] * (M_PI / 180.0);
		double Lrad = L[i] * (M_PI / 180.0);
		X[i] = (R + H[i]) * cos(Brad) * cos(Lrad);
		Y[i] = (R + H[i]) * cos(Brad) * sin(Lrad);
		Z[i] = (R + H[i]) * sin(Brad);

		x[i] = (R + H[i] - D) * cos(Brad) * cos(Lrad);
		y[i] = (R + H[i] - D) * cos(Brad) * sin(Lrad);
		z[i] = (R + H[i] - D) * sin(Brad);

		nx[i] = cos(B[i] * (M_PI / 180.0)) * cos(Lrad);
		ny[i] = cos(B[i] * (M_PI / 180.0)) * sin(Lrad);
		nz[i] = sin(B[i] * (M_PI / 180.0));
	}

# pragma omp parallel for private(j,rx,ry,rz,r_ij)
	for ( i = 0; i < N; i++)
	{
		for ( j = 0; j < N; j++)
		{
			rx = X[i] - x[j];
			ry = Y[i] - y[j];
			rz = Z[i] - z[j];
			r_ij = sqrt(rx * rx + ry * ry + rz * rz);


			A[i][j] = (1.0 / (4.0 * M_PI * r_ij * r_ij * r_ij)) * (rx * nx[i] + ry * ny[i] + rz * nz[i]);
		}
	}


	//for (int i = 0; i < 5; i++)
	//{
	//    printf("%d: %.20lf\n",myrank,Alocal[i][i]);
	//}

	// BiCG Stab

	xj = new double[N]();
	r0 = new double[N]();
	r_prev = new double[N]();
	pj = new double[N]();
	x_prev = new double[N]();
	rozdiel = new double[N]();
	s = new double[N]();
	rj = new double[N]();

	for ( i = 0; i < N; i++) {
		xj[i] = 0.0;
		r0[i] = g[i];
		r_prev[i] = g[i];
		pj[i] = r_prev[i];
		s[i] = 0.0;
		rj[i] = 0.0;
	}

	temp_result = new double[N]();
	temp_result2 = new double[N]();

	double citatel, menovatel, norma;

	do {
		for ( i = 0; i < N; i++) {
			x_prev[i] = xj[i];
		}

# pragma omp parallel for private(j,sum)
		for ( i = 0; i < N; i++) {
			sum = 0.0;
			for ( j = 0; j < N; j++) {
				sum += A[i][j] * pj[j];
			}
			Apj[i] = sum;
		}

		citatel = 0, menovatel = 0;
		for ( i = 0; i < N; i++) {
			citatel += r_prev[i] * r0[i];
		}
		for ( i = 0; i < N; i++) {
			menovatel += Apj[i] * r0[i];
		}
		alfa = citatel / menovatel;

		for ( i = 0; i < N; i++) {
			temp_result[i] = alfa * Apj[i];
		}
		for ( i = 0; i < N; i++) {
			s[i] = r_prev[i] - temp_result[i];
		}

# pragma omp parallel for private(j,sum)
		for ( i = 0; i < N; i++) {
			sum = 0.0;
			for ( j = 0; j < N; j++) {
				sum += A[i][j] * s[j];
			}
			As[i] = sum;
		}

		citatel = 0, menovatel = 0;
		for ( i = 0; i < N; i++) {
			citatel += As[i] * s[i];
		}
		for ( i = 0; i < N; i++) {
			menovatel += As[i] * As[i];
		}
		omega = citatel / menovatel;

		for ( i = 0; i < N; i++) {
			temp_result[i] = alfa * pj[i];
		}
		for ( i = 0; i < N; i++) {
			temp_result[i] = xj[i] + temp_result[i];
		}
		for ( i = 0; i < N; i++) {
			temp_result2[i] = omega * s[i];
		}
		for ( i = 0; i < N; i++) {
			xj[i] = temp_result[i] + temp_result2[i];
		}


		for ( i = 0; i < N; i++) {
			temp_result[i] = omega * As[i];
		}
		for ( i = 0; i < N; i++) {
			rj[i] = s[i] - temp_result[i];
		}


		citatel = 0, menovatel = 0;
		for (int i = 0; i < N; i++) {
			citatel += rj[i] * r0[i];
		}
		for ( i = 0; i < N; i++) {
			menovatel += r_prev[i] * r0[i];
		}
		beta = (citatel / menovatel) * (alfa / omega);


		for ( i = 0; i < N; i++) {
			temp_result[i] = omega * Apj[i];
		}
		for ( i = 0; i < N; i++) {
			temp_result[i] = pj[i] - temp_result[i];
		}
		for ( i = 0; i < N; i++) {
			temp_result[i] = beta * temp_result[i];
		}
		for ( i = 0; i < N; i++) {
			pj[i] = rj[i] + temp_result[i];
		}


		for ( i = 0; i < N; i++) {
			r_prev[i] = rj[i];
		}

		norma = 0.0;

		for ( i = 0; i < N; i++) {
			norma += rj[i] * rj[i];
		}
		norma = sqrt(norma);

		iter++;
		std::cout << "iter: " << iter << " " << norma << std::endl;

	} while (iter < maxIterations && norma > epsilon);

	double G;
# pragma omp parallel for private(j,G,rx,ry,rz,r_ij)
	for ( i = 0; i < N; i++)
	{
		G = 0;
		for ( j = 0; j < N; j++)
		{
			rx = X[i] - x[j];
			ry = Y[i] - y[j];
			rz = Z[i] - z[j];
			r_ij = sqrt(rx * rx + ry * ry + rz * rz);
			G += 1.0 / (4.0 * M_PI * r_ij) * xj[j];
		}
		u[i] = G;
	}

	file = fopen("output_omp.dat", "w");
	if (!file) {
		perror("Error opening file");
	}

	for ( i = 0; i < N; i++) {
		fprintf(file, "%lf %lf %lf\n", B[i], L[i], u[i]);
	}

	fclose(file);

	delete[] B;
	delete[] L;
	delete[] H;
	delete[] g;
	delete[] u2n2;
	delete[] Apj;
	delete[] As;
	delete[] X;
	delete[] Y;
	delete[] Z;
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] nx;
	delete[] ny;
	delete[] nz;

	delete[] xj;
	delete[] r0;
	delete[] r_prev;
	delete[] pj;
	delete[] x_prev;
	delete[] rozdiel;
	delete[] s;
	delete[] rj;
	delete[] temp_result;
	delete[] temp_result2;

	for (int i = 0; i < N; i++) {
		delete[] A[i];
	}
	delete[] A;

	delete[] u;

}