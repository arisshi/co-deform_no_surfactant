/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <math.h>
#include <string>
#include "../include/fd.h"
#include "../include/rw.h"

using namespace std;

int main(){
	int i, j, m, n;
	char fn[1024];
	double params[50] = {0.0};

	// read parameters
	sprintf(fn, "./params.in" );
	readInput(fn, params);

	double Ca      = params[0];
	double Fr      = params[1];
	double Ma      = params[2];
	double tstop   = params[3];
	double r1      = params[4];
	double t1      = params[5];
	double dr      = params[6];
	double dt      = params[7];
	double dtrec   = params[8];
	double po      = params[9];
	double pobar   = params[10];
	
	int J  = r1/dr;
	int N  = t1/dt + 1;
	int M = N*dt/dtrec;
	
	int J1 = J+1;
	int M1 = M+1;
	int M1J1 = M1*J1;
	int M1J  = M1*J;
	
	// allocate memory and initialize to zero
	double *R  = (double*) calloc(M1J1, sizeof(double));
	double *T  = (double*) calloc(M1J1, sizeof(double));
	double *H  = (double*) calloc(M1J1, sizeof(double));
	double *E  = (double*) calloc(M1J1, sizeof(double));
	double *F  = (double*) calloc(M1J1, sizeof(double));
	double *B  = (double*) calloc(M1J1, sizeof(double));
	double *Q  = (double*) calloc(M1J1, sizeof(double));
	double *P  = (double*) calloc(M1J1, sizeof(double));
	double *gamma1  = (double*) calloc(M1J1, sizeof(double));
	double *gamma2  = (double*) calloc(M1J1, sizeof(double));
	double *h1    = (double*) calloc(M1J1, sizeof(double));
	double *h2    = (double*) calloc(M1J1, sizeof(double));
	double *pdyn    = (double*) calloc(M1J1, sizeof(double));
	double *v1      = (double*) calloc(M1J1, sizeof(double));
	double *v2      = (double*) calloc(M1J1, sizeof(double));

	//printf("before fdevol IN RUN.CPP.");
	// time evolution of h, g, and f
	fdevol(J, N, M, dr, dt, params, H, E, F, Q, P);

	//printf("before fdaux IN RUN.CPP.");
	// back-calculate p, q, and vs
     //fdaux(J, N, M, dr, dt, params, H, E, F, Q, gamma1, gamma2, h1, h2, pdyn, v1, v2);

	//printf("before fdgrid IN RUN.CPP.");
	// space and time domains
	fdgrid(J, N, M, dr, dt, R, T);

	for (i = 0; i < M1J1; i++){B[i] = F[i] + E[i];}

	//printf("THIS IS BEFORE THE SPRINTF IN RUN.CPP.");
	// write to file
	sprintf(fn, "evol_r" );      write(J, M, dr, dt, params, R , fn);
	sprintf(fn, "evol_t" );      write(J, M, dr, dt, params, T , fn);
	sprintf(fn, "evol_h" );      write(J, M, dr, dt, params, H , fn);
	sprintf(fn, "evol_e" );      write(J, M, dr, dt, params, E , fn);
	sprintf(fn, "evol_f" );      write(J, M, dr, dt, params, F , fn);
	sprintf(fn, "evol_b" );      write(J, M, dr, dt, params, B , fn);
	sprintf(fn, "evol_p" );      write(J, M, dr, dt, params, P , fn);
	sprintf(fn, "evol_q" );      write(J, M, dr, dt, params, Q , fn);
//	sprintf(fn, "aux_gamma1" ); write(J, M, dr, dt, params, gamma1 , fn);
//	sprintf(fn, "aux_gamma2" ); write(J, M, dr, dt, params, gamma2 , fn);
//	sprintf(fn, "aux_h1" );   write(J, M, dr, dt, params, h1 , fn);
//	sprintf(fn, "aux_h2" );   write(J, M, dr, dt, params, h2 , fn);
//	sprintf(fn, "aux_pdyn" );   write(J, M, dr, dt, params, pdyn , fn);
//	sprintf(fn, "aux_v1");      write(J, M, dr, dt, params, v1, fn);
//	sprintf(fn, "aux_v2");      write(J, M, dr, dt, params, v2, fn);

	//printf("AFTER THE SPRINTF IN RUN.CPP.");

	// free memory
	free(R);
	free(T);
	free(H);
	free(E);
	free(F);
	free(B);
	free(Q);
	free(P);
	free(gamma1);
	free(gamma2);
	free(h1  );
	free(h2  );
	free(pdyn  );
	free(v1    );
	free(v2    );
	
	return(0);
}
