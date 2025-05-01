#define _CRT_SECURE_NO_WARNINGS
#include <mpi.h>
#include <iostream>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define N 160002
#define R 6378000.0
#define GM 398600.5
#define D 50000//o kolko je vnoreny bod


int main(int argc, char** argv) {

    int nprocs, myrank, nlocal, istart, iend, nlast, iglobal;
    double r_ij, rx, ry, rz, skalarny;

    double* B = new double[N]();
    double* L = new double[N]();
    double* H = new double[N]();
    double* g = new double[N]();
    double* u2n2 = new double[N]();
    double* temp_g;
    double G;
    double** Alocal;
    double* Apj = new double[N + 100]();
    double* As = new double[N + 100]();
    double* u_temp;


    double* Apj_temp;
    double* As_temp;

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

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    if (myrank == 0)    //citanie suboru v procese 0
    {
        FILE* file = fopen("BL-160002.dat", "r");

        if (!file) {
            perror("Error opening file");
        }

        for (int i = 0; i < N; i++)
        {
            int count = fscanf(file, "%lf %lf %lf %lf %lf", &B[i], &L[i], &H[i], &g[i], &u2n2[i]);
            //g[i] = GM / (R * R);
            //H[i] = 0;
            g[i] = g[i] * 0.00001;


        }
        fclose(file);
    }

    nlocal = (int)(N / nprocs + 1);

    temp_g = new double[N]();
    Alocal = new double* [nlocal];
    nlast = N - nlocal * (nprocs - 1);

    for (int i = 0; i < nlocal; i++)
    {
        Alocal[i] = new double[N]();
    }

    //rozposlanie vektorov vsetkym procesom
    MPI_Bcast(B, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(L, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(H, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(u2n2, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    istart = 0;

    if (myrank == (nprocs - 1))
    {
        iend = nlast;
    }
    else
    {
        iend =  nlocal;
    }

    //printf("\n%d: %d-%d\n\n", myrank, istart, iend);

    for (int i = 0; i < N; i++)
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


    for (int i = 0; i < iend; i++)
    {
        for (int j = 0; j < N; j++)
        { 
            iglobal = myrank * nlocal + i;
            rx = X[iglobal] - x[j];
            ry = Y[iglobal] - y[j];
            rz = Z[iglobal] - z[j];
            r_ij = sqrt(rx * rx + ry * ry + rz * rz);

            //skalarny sucin
            skalarny = rx * nx[iglobal] + ry * ny[iglobal] + rz * nz[iglobal];

            Alocal[i][j] = (1.0 / (4.0 * M_PI * r_ij * r_ij * r_ij)) * skalarny;
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

    for (int i = 0; i < N; i++) {
        xj[i] = 0.0;
        r0[i] = g[i];
        r_prev[i] = g[i];
        pj[i] = r_prev[i];
        s[i] = 0.0;
        rj[i] = 0.0;
    }

    temp_result = new double[N]();
    temp_result2 = new double[N]();
    Apj_temp = new double[N]();
    As_temp = new double[N]();


    double citatel, menovatel, norma;

    do {
        for (int i = 0; i < N; i++) {
            x_prev[i] = xj[i];
        }

        for (int i = 0; i < iend; i++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++) {
                sum += Alocal[i][j] * pj[j];
            }
            Apj_temp[i] = sum;
        }
        MPI_Allgather(Apj_temp, nlocal, MPI_DOUBLE, Apj, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

        citatel = 0, menovatel = 0;
        for (int i = 0; i < N; i++) {
            citatel += r_prev[i] * r0[i];
        }
        for (int i = 0; i < N; i++) {
            menovatel += Apj[i] * r0[i];
        }
        alfa = citatel / menovatel;

        for (int i = 0; i < N; i++) {
            s[i] = r_prev[i] - alfa * Apj[i];
        }


        for (int i = 0; i < iend; i++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++) {
                sum += Alocal[i][j] * s[j];
            }
            As_temp[i] = sum;
        }
        MPI_Allgather(As_temp, nlocal, MPI_DOUBLE, As, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

        citatel = 0, menovatel = 0;
        for (int i = 0; i < N; i++) {
            citatel += As[i] * s[i];
            menovatel += As[i] * As[i];
        }
        omega = citatel / menovatel;

        for (int i = 0; i < N; i++) {
            xj[i] = xj[i] + alfa * pj[i] + omega * s[i];
        }

        for (int i = 0; i < N; i++) {
            rj[i] = s[i] - omega * As[i];
        }


        citatel = 0, menovatel = 0;
        for (int i = 0; i < N; i++) {
            citatel += rj[i] * r0[i];
        }
        for (int i = 0; i < N; i++) {
            menovatel += r_prev[i] * r0[i];
        }
        beta = (citatel / menovatel) * (alfa / omega);


        for (int i = 0; i < N; i++) {
            pj[i] = rj[i] + beta * (pj[i] - omega * Apj[i]);
        }


        for (int i = 0; i < N; i++) {
            r_prev[i] = rj[i];
        }

        norma = 0.0;

        for (int i = 0; i < N; i++) {
            norma += rj[i] * rj[i];
        }
        norma = sqrt(norma);

        iter++;
        if (myrank == 0)std::cout << myrank << ": " << "iter: " << iter << " " << norma << std::endl;

    } while (iter < maxIterations && norma > epsilon);

    if (myrank == 0) {
        FILE* file = fopen("alfy.txt", "w");

        if (!file) {
            perror("Error opening file");
        }

        for (int i = 0; i < N; i++) {
            fprintf(file, "%lf\n", xj[i]);
        }

        fclose(file);
    }

    u_temp = new double[nlocal]();

    for (int i = 0; i < iend; i++)
    {
        G = 0;
        for (int j = 0; j < N; j++)
        {
            iglobal = myrank * nlocal + i;
            rx = X[iglobal] - x[j];
            ry = Y[iglobal] - y[j];
            rz = Z[iglobal] - z[j];
            r_ij = sqrt(rx * rx + ry * ry + rz * rz);
            G += 1.0 / (4.0 * M_PI * r_ij) * xj[j];
        }
        u_temp[i] = G;
    }
    MPI_Gather(u_temp, nlocal, MPI_DOUBLE, u, nlocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (myrank == 0) {
        FILE* file = fopen("output.dat", "w");
        if (!file) {
            perror("Error opening file");
        }

        for (int i = 0; i < N; i++) {
            fprintf(file, "%lf %lf %lf\n", B[i], L[i], u[i]);
        }

        fclose(file);
    }
    MPI_Finalize();

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

    delete[] Apj_temp;
    delete[] As_temp;
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

    for (int i = 0; i < nlocal; i++) {
        delete[] Alocal[i];
    }
    delete[] Alocal;
    delete[] temp_g;
    
    delete[] u_temp;

    //delete[] u; //????????????
}