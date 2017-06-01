/* FINITE DIFFERENCE METHOD
 *  Finite difference method for the one-dimensional evolution of film depth and insoluble surfactant
 *  concentration on an initially planar free surface penetrated by a rigid sphere from below.
 *
 *  The space and time domains are discretized on a uniform stencil with grid spacing dr and dt,
 *  respectively. Cylindrical polar coordinates are used in space. Time advancement is carried out
 *  using a trapezoidal method with centered-difference spatial operators (Crank-Nicholson scheme).
 *
 * REFERENCES
 *  von Rosenberg, Methods for the Numerical Solution of Partial Differential Equations (Elsevier, 1969)
 *  
 * PARAMETERS
 *  J		[input]			number of space segments
 *  N		[input]			number of time segments calculated
 *  M		[input]			number of time segments recorded
 *  dr  [input]			spatial grid spacing
 *  dt  [input]			temporal grid spacing
 *
 *  h   [output]		film thickness
 *  g   [output]		excess surface concentration
 *  f   [output]		sphere position
 *  p   [output]		dynamic pressure
 *  q   [output]		volume flux per unit length
 *  vs  [output]		surface velocity
 */

#ifndef FD_H
#define FD_H


/* HEADER FILES */
#include "la.h"
#include "bessel.h"
#include <Eigen/Sparse>
#include <vector>
#include <math.h>


using namespace std;

/* PROTOTYPES */
void fdgrid(int, int, int, double, double, double *, double *);
void fdmat (int, double, double *, double *, double *, double *, double *, double *, double *, double *, double *);
void fdevol(int, int, int, double, double, double *, double *, double *, double *, double *, double *);
void fdstep(int, int, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
//void fddiff(int, int, double , double, double,
//            double *, double *, double *, 
//            double *, double *, double * );
//void fddiffhalf(int , int , double, double, double,
//            double *, double *, double *,
//            double *, double *, double *, double *);
//void fdaux(int, int, int, double, double, double *,
//           double *, double *, double *, double *,
//		 double *, double *, double *, double *, double *, double *, double *);



/* IMPLEMENTATIONS */
// space and time stencils
void fdgrid(int J, int N, int M, double dr, double dt, double *R, double *T){
	int i, n, m;
	int J1 = J+1;
	int N1 = N+1;
	int M1 = M+1;
	int nrec = N/M;
	double r, t;

	m = 0;
	for (n = 0; n < N1; n++){
		t = n*dt;
		if (n % nrec == 0){
			for (i = 0; i < J1; i++){
				r = i*dr;
				R[m*J1 + i] = r;
				T[m*J1 + i] = t;
			}
			m++;
		}
	}
}

// unit and derivative matrices
void fdmat(int J, double dr,  double *U, double *DF, double *DB, double *DC, double *DA, double *DD, double *DP, double *DN, double *DM){

	int i, j;
	int JJ = J*J;
	double dr2 = dr*dr;
	double dr4 = dr2*dr2;

	// initialize
	for (i = 0; i < JJ; i++){
		U [i] = 0.0;
		DF[i] = 0.0;
		DB[i] = 0.0;
		DC[i] = 0.0;
		DA[i] = 0.0;
		DD[i] = 0.0;
		DP[i] = 0.0;
		DN[i] = 0.0;
		DM[i] = 0.0;
	}

	// identity matrix
	for (i = 0; i < J; i++){
		U[i*J + i] = 1;
	}

	// derivative matrice1
	for (i = 1; i < J; i++){
		// forward
		DF[ i   *J +  i   ] = -1;
		DF[(i-1)*J +  i   ] =  1;

		// backward
		DB[ i   *J +  i   ] =  1;
		DB[ i   *J + (i-1)] = -1;

		// centered
		DC[ i   *J + (i-1)] = -0.5;
		DC[(i-1)*J +  i   ] =  0.5;

		// backward-gradient
		DA[ i   *J +  i   ] =  1 + 0.5/i;
		DA[ i   *J + (i-1)] = -1 + 0.5/i;

		// positive convection speed
		DP[i*J+i] = 1;
		DP[i*J+i-1] = 1/i-1;
	
		// negative convection speed
		DN[i*J+i] = -1;
		if (i < J-1) {DN[i*J+i+1] = 1+1/i;}
	
	
		DM[i*J+i] = 0.5;
		if (i < J-1) {DM[i*J+i+1] = 0.5;}
	}
	DF[0] = -1; DF[1] =  1;
	DB[0] =  1; DB[1] = -1;
	DC[1] =  0; DC[(J-1)*J + (J-1)] =  1.5; DC[(J-1)*J + (J-2)] = -2; DC[(J-1)*J + (J-3)] = 0.5;
	DA[0] =  4; 
	DP[0] =  2; DP[1] = -2; 
	DN[0] = -2; DN[1] =  2;
	DM[0] = 0.5; DM[1] = 0.5;
		
	// centered-Laplacian
	// METHOD #1: input entries manually
//	for (i = 1; i < J-1; i++){
//		DD[ i   *J +  i   ] = -2;
//		DD[ i   *J + (i-1)] =  1 - 0.5/i;
//		DD[ i   *J + (i+1)] =  1 + 0.5/i;
//	}
//	DD[0] = -4; DD[1] = 4; DD[(J-1)*J + (J-1)] = -2; DD[(J-1)*J + (J-2)] = 1 - 0.5/(J-1);

//
//	// METHOD #2: matrix multiplication
//	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, J, J, J, 1.0, DA, J, DF, J, 0.0, DD, J);
	// METHOD #3: banded matrix multiplication
	laband(J, DA, 0, 1, DF, 1, 0, DD);	


	


	for (i = 0; i < JJ; i++){
		DF[i] /= dr;
		DB[i] /= dr;
		DC[i] /= dr;
		DA[i] /= dr;
		DD[i] /= dr2;
		DP[i] /= dr;
		DN[i] /= dr;
	}
//	for (i = 0; i < J*J; i++){
//		printf("DD(%i) = %f; \n", i+1, DD[i]);
//	}

}






// advance +1 timestep
//  h0, h1 are (J+1)-vectors
void fdstep(int J, int n, double *params, double *h0, double *e0, double *f0, double *q0, double *h1, double *e1, double *f1, double *q1, double *p1){

typedef Eigen::SparseMatrix<double> spmat;
typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(18*J);

	int i, j;
	int JJ = J*J;
	int J2 = J + J;
	int J3 = J2+ J;

	double Ca    = params[0];
	double Bo    = params[1];
	double Ma    = params[2];
	double tstop = params[3];
	double R1    = params[4];
	double dr    = params[6];
	double dt    = params[7];

	double dr2  = dr*dr;
	double dr4  = dr2*dr2;
	double time = tstop; if (n*dt <= tstop) {time = n*dt;}

	// allocate memory and initialize to zero

	double *A    = (double*) calloc(JJ, sizeof(double));
	double *B    = (double*) calloc(JJ, sizeof(double));
	double *C    = (double*) calloc(JJ, sizeof(double));
	double *D    = (double*) calloc(JJ, sizeof(double));
	double *E    = (double*) calloc(JJ, sizeof(double));
	double *DA    = (double*) calloc(JJ, sizeof(double));
	double *DD    = (double*) calloc(JJ, sizeof(double));
	double *DF    = (double*) calloc(JJ, sizeof(double));
	double *DB    = (double*) calloc(JJ, sizeof(double));
	double *DC    = (double*) calloc(JJ, sizeof(double));
	double *DN    = (double*) calloc(JJ, sizeof(double));
	double *DM    = (double*) calloc(JJ, sizeof(double));
	double *DP    = (double*) calloc(JJ, sizeof(double));
	double *U     = (double*) calloc(JJ, sizeof(double));
	
	double *L11  = (double*) calloc(JJ, sizeof(double));
	double *L12  = (double*) calloc(JJ, sizeof(double));
	double *L13 = (double*) calloc(JJ, sizeof(double));
	double *L21  = (double*) calloc(JJ, sizeof(double));
	double *L22  = (double*) calloc(JJ, sizeof(double));
	double *L23  = (double*) calloc(JJ, sizeof(double));
	double *L22p = (double*) calloc(JJ, sizeof(double));
	double *L33  = (double*) calloc(JJ, sizeof(double));
	double *L34  = (double*) calloc(JJ, sizeof(double));
	double *L31  = (double*) calloc(JJ, sizeof(double));
	double *L32  = (double*) calloc(JJ, sizeof(double));
	double *zero = (double*) calloc(JJ, sizeof(double));

	double *b1   = (double*) calloc(J , sizeof(double));
	double *b2   = (double*) calloc(J , sizeof(double));
	double *b3   = (double*) calloc(J , sizeof(double));
	double *b4   = (double*) calloc(J , sizeof(double));
	double *ones = (double*) calloc(J , sizeof(double));
	double *anh  = (double*) calloc(J , sizeof(double));
	double *tp1  = (double*) calloc(J , sizeof(double));
	double *tp2  = (double*) calloc(J , sizeof(double));
	double *tp3  = (double*) calloc(J , sizeof(double));
	double *tp4  = (double*) calloc(J , sizeof(double));
	double *tp5  = (double*) calloc(J , sizeof(double));
	
	double *hs    = (double*) calloc(J , sizeof(double));
	double *es    = (double*) calloc(J , sizeof(double));
	double *fm    = (double*) calloc(J , sizeof(double));
	double *qs    = (double*) calloc(J , sizeof(double));
	double *hm   = (double*) calloc(J , sizeof(double));
	double *em   = (double*) calloc(J , sizeof(double));
	double *hh   = (double*) calloc(J , sizeof(double));
	double *qm   = (double*) calloc(J , sizeof(double));
	double *pm = (double*) calloc(J , sizeof(double));
	double *thm = (double*) calloc(J , sizeof(double));
	double *h20m = (double*) calloc(J , sizeof(double));
	double *vr = (double*) calloc(J , sizeof(double));
     int ncut =-1;

	fdmat(J, dr, U, DF, DB, DC, DA, DD, DP, DN, DM);
	
	// h is h1, e is Ca*h21, f is h20, q is radial veclocity, p is pressure
	for (i = 0; i < J+1; i++){f1[i] = -1.0 - i*i*dr2/2.0 + time;}

	for (i = 0; i < J; i++){	
		hm[i] = (h0[i]+h0[i+1])/2.0;	
		qm[i] = (q0[i]+q0[i+1])/2.0;
		em[i] = (e0[i]+e0[i+1])/2.0;
		fm[i] = (f0[i]+f0[i+1])/2.0;
		thm[i]  = hm[i] - em[i] - fm[i];
		h20m[i] = fm[i] + dt;
	}


	// assemble block matrices
	for (i = 0; i < J; i++){	A[i*J+i] = qm[i];}
	laband(J,DA,1,1,A ,1,1,B);
	laband(J,B ,1,1,DM,1,1,C);

	for (i = 0; i < JJ; i++){
		L11[i] =  C[i]*dt/2.0;
	if(n > ncut){L12[i] = -C[i]*dt/2.0;}
		L21[i] = DD[i];
		L23[i] = 3.0*Ca*DA[i];
		L32[i] = DD[i];
		L33[i] = -3.0*Ca*DA[i];
	}
	for (i = 0; i < J; i++){
		L11[i*J+i] += 1.0;
	if(n>ncut){ L12[i*J+i] +=-1.0;}
		L21[i*J+i] -= Bo;
		L32[i*J+i] += Bo;
	}
	// Robin condition at R1
	i = J - 1;
	double m = sqrt(Bo);
	double gridfactor = (1.0+1.0/2.0/i);
	double myy = m*bessy1(m*R1)/bessy0(m*R1);
	L32[i*J+i  ] += 1.0/dr2*gridfactor*(-2.0*dr)*myy;
	L32[i*J+i-1] += 1.0/dr2*gridfactor;
	if (n>ncut){
	     L12[i*J+i  ] += dt/2.0/dr*gridfactor*qm[i]*(dr*myy);
		L12[i*J+i-1] += dt/2.0/dr*gridfactor*qm[i]*(-1.0/2.0);
	}

	// assemble right hand side
	// Freidrichs scheme
	laband(J,B,1,1,thm,tp1);
	laband(J,B,1,1,h20m,tp2);
	if (n <= ncut){
		for(i = 0; i < J; i++){b1[i] = dt - dt/2.0*(tp1[i]-tp2[i]) + h0[i];}
	} else if (n > ncut and n*dt <= tstop) {
		for(i = 0; i < J; i++){b1[i] = dt - dt/2.0*(tp1[i]-tp2[i]) + (h0[i]-e0[i]);}
	}else{
		for(i = 0; i < J; i++){b1[i] =    - dt/2.0*(tp1[i]-tp2[i]) + (h0[i]-e0[i]);}
	}
	for (i = 0; i < J; i++){ b2[i] = 0.0; b3[i] = 0.0; }


	// put the matrices into the sparse system
	for(i = 0; i < J; i++){
		if(i == 0){
			tripletList.push_back( T(i+J , i     , L21[i*J+i  ]));
			tripletList.push_back( T(i+J , i+1   , L21[i*J+i+1]));
			
			tripletList.push_back( T(i+J , i+J2  , L23[i*J+i  ]));
			tripletList.push_back( T(i+J , i+J2+1, L23[i*J+i+1]));
			
			tripletList.push_back( T(i+J2, i+J   , L32[i*J+i  ]));
			tripletList.push_back( T(i+J2, i+J+1 , L32[i*J+i+1]));
			
			tripletList.push_back( T(i+J2, i+J2  , L33[i*J+i  ]));
			tripletList.push_back( T(i+J2, i+J2+1, L33[i*J+i+1]));
		
			tripletList.push_back( T(i   , i     , L11[i*J+i])   );
			tripletList.push_back( T(i   , i+1   , L11[i*J+i+1]));

			tripletList.push_back( T(i   , i+J   , L12[i*J+i])   );
			tripletList.push_back( T(i   , i+J+1 , L12[i*J+i+1]));
		}else if(i == J-1){
			tripletList.push_back( T(i+J , i     , L21[i*J+i  ]));
			tripletList.push_back( T(i+J , i-1   , L21[i*J+i-1]));
			
			tripletList.push_back( T(i+J , i+J2  , L23[i*J+i  ]));
			tripletList.push_back( T(i+J , i+J2-1, L23[i*J+i-1]));
			
			tripletList.push_back( T(i+J2, i+J   , L32[i*J+i  ]));
			tripletList.push_back( T(i+J2, i+J-1 , L32[i*J+i-1]));
			
			tripletList.push_back( T(i+J2, i+J2  , L33[i*J+i  ]));
			tripletList.push_back( T(i+J2, i+J2-1, L33[i*J+i-1]));
		
			tripletList.push_back( T(i   , i     , L11[i*J+i])   );
			tripletList.push_back( T(i   , i-1   , L11[i*J+i-1]));

			tripletList.push_back( T(i   , i+J   , L12[i*J+i])   );
			tripletList.push_back( T(i   , i+J-1 , L12[i*J+i-1]));
		}else{
			tripletList.push_back( T(i+J , i     , L21[i*J+i  ]));
			tripletList.push_back( T(i+J , i+1   , L21[i*J+i+1]));
			tripletList.push_back( T(i+J , i-1   , L21[i*J+i-1]));
			
			tripletList.push_back( T(i+J , i+J2  , L23[i*J+i  ]));
			tripletList.push_back( T(i+J , i+J2+1, L23[i*J+i+1]));
			tripletList.push_back( T(i+J , i+J2-1, L23[i*J+i-1]));
			
			tripletList.push_back( T(i+J2, i+J   , L32[i*J+i  ]));
			tripletList.push_back( T(i+J2, i+J+1 , L32[i*J+i+1]));
			tripletList.push_back( T(i+J2, i+J-1 , L32[i*J+i-1]));
			
			tripletList.push_back( T(i+J2, i+J2  , L33[i*J+i  ]));
			tripletList.push_back( T(i+J2, i+J2+1, L33[i*J+i+1]));
			tripletList.push_back( T(i+J2, i+J2-1, L33[i*J+i-1]));
		
			tripletList.push_back( T(i   , i     , L11[i*J+i])   );
			tripletList.push_back( T(i   , i+1   , L11[i*J+i+1]));
			tripletList.push_back( T(i   , i-1   , L11[i*J+i-1]));

			tripletList.push_back( T(i   , i+J   , L12[i*J+i])   );
			tripletList.push_back( T(i   , i+J+1 , L12[i*J+i+1]));
			tripletList.push_back( T(i   , i+J-1 , L12[i*J+i-1]));
		}
	}






	Eigen::SparseMatrix<double> Lall(J3, J3);
	Lall.setFromTriplets(tripletList.begin(),tripletList.end());



	// polulate the right hand side 
	Eigen::VectorXd ball(J3);
	
	for(i = 0; i < J; i++)   {ball(i) = b1[i];}
	for(i = J; i < J2; i++)   {ball(i) = b2[i-J];}
	for(i = J2; i < J3; i++)   {ball(i) = b3[i-J2];}


	// solve linear system
	//Eigen::SimplicialCholesky<spmat> solver(Lall);
	Eigen::SparseLU<spmat, Eigen::COLAMDOrdering<int> > solver(Lall);
	//Eigen::SparseLU<spmat, Eigen::AMDOrdering<int> > solver(Lall);
	Eigen::VectorXd x = solver.solve(ball);
	
	// copy solution and assign right BCs
	for(i = 0; i < J; i++)	     h1[i]      = x(i);
	for(i = J; i < J2; i++)	     e1[i-J]    = x(i);
	h1[J] = 0.0;
	e1[J] = e1[J-2] - 2.0*dr*myy*e1[J-1];
	
	
	for(i = J2; i < J3-1; i++)	q1[i-J2+1] = (x(i)+x(i+1))/2.0;
	q1[0] = 0.0;
	q1[J] = x(J3-1);
	
	//q1[0] = 0.0;
	//for(i = J2; i < J3-1; i++)	q1[i-J2+1] = x(i);
	
	

	// compute pressure 
	for(i = J2; i < J3; i++)	tp1[i-J2] = x(i);
	laband(J,DA,0,1,tp1,pm);
	for(i = 0; i < J-1; i++)	p1[i+1] = (pm[i] + pm[i+1])/2.0;
	p1[0] = pm[0];
	p1[J] = 0.0;

//	// debug print Lall * x = ball
//	for (int k = 0; k < Lall.outerSize(); ++k){
//		for(Eigen::SparseMatrix<double>::InnerIterator it(Lall,k); it; ++it){
//			std::cout << "Lall(" << it.row()+1 << ", ";
//			std::cout << it.col()+1 << " ) = " << it.value() << "; \n";
//		}
//	}
//	for (i = 0; i < J3; i++){
//		printf("ball(%i) = %f; \n", i+1, ball(i));
//	}
//	for (i = 0; i < J3; i++){
//		printf("x(%i) = %f; \n", i+1, x(i));
//	}
	
	

	// free memory
	free(A    );
	free(B    );
	free(C    );
	free(D    );
	free(E    );
	free(DA   );
	free(DD   );
	free(L11  );
	free(L13 );
	free(L12  );
	free(L21  );
	free(L22  );
	free(L23  );
	free(L22p );
	free(L33  );
	free(L34  );
	free(L31  );
	free(L32  );
	free(zero );
	free(b1   );
	free(b2   );
	free(b3   );
	free(b4   );
	free(anh  );
	free(ones );
	free(tp1  );
	free(tp2  );
	free(tp3  );
	free(tp4  );
	free(tp5  );
	free(hs    );
	free(es    );
	free(fm    );
	free(qs    );
	free(hm   );
	free(em   );
	free(hh   );
	free(qm   );
	free(pm);
	free(thm);
	free(h20m);
	free(vr);
	free(DF);
	free(DB);
	free(DC);
	free(DN);
	free(DM);
	free(DP);
	free(U );

}


// time evolution
void fdevol(int J, int N, int M, double dr, double dt, double *params,
					  double *H, double *E, double *F, double *Q, double *P){
	int i, j, n, m, nrec;
	int J1 = J+1;
	int N1 = N+1;
	int M1 = M+1;
	double t, r;
	double dr2 = dr*dr;

	double Ca = params[0];
	double Bo = params[1];
	
	// allocate memory and initialize to zero
	double *h0 = (double*) calloc(J1, sizeof(double));
	double *e0 = (double*) calloc(J1, sizeof(double));
	double *f0 = (double*) calloc(J1, sizeof(double));
	double *q0 = (double*) calloc(J1, sizeof(double));
	double *h1 = (double*) calloc(J1, sizeof(double));
	double *e1 = (double*) calloc(J1, sizeof(double));
	double *f1 = (double*) calloc(J1, sizeof(double));
	double *q1 = (double*) calloc(J1, sizeof(double));
	double *p0 = (double*) calloc(J1, sizeof(double));
	double *p1 = (double*) calloc(J1, sizeof(double));

	double filmsum = 0.0;
	double temp;

	// initialize record counter and wait time
	m = 0;
	nrec = N/M;
	// initial conditions
	for (i = 0; i < J1; i++){

		H[0*J1+ i] = 0.0;
		E[0*J1+ i] = 0.0;
		Q[0*J1+ i] = 0.0; // Initial guess for q at t=0.
		P[0*J1+ i] = 0.0; 
		F[0*J1+ i] = -1.0-dr2*i*i/2.0; 
		
		h0[i] = H[0*J1+ i];
		q0[i] = Q[0*J1+ i];
		e0[i] = E[0*J1+ i];
		f0[i] = F[0*J1+ i];

	}
	

	// time-evolve
	m = 0;
	for (n = 1; n < N1; n++){

		fdstep(J, n, params, h0, e0,f0, q0, h1, e1,f1, q1, p1);

//		filmsum = 0.0;
		for (i = 0; i < J1; i++){
//			temp = h1[i] - (f1[i] + e1[i]);
//			if (temp < 0.0){filmsum += temp;}
			h0[i] = h1[i];
			q0[i] = q1[i];
			e0[i] = e1[i];
			f0[i] = f1[i];
		}

//		if (filmsum < 0.0){
//			printf("finalrow = %i, n = %i / %i Simulation stops here \n", m, n, N);
//			m++;
//			for (i = 0; i < J1; i++){
//				H[m*J1+ i] = h1[i];
//				Q[m*J1+ i] = q1[i];
//				E[m*J1+ i] = e1[i];
//				P[m*J1+ i] = p1[i];
//				F[m*J1+ i] = f1[i];
//			}
//			return;
//		}

		if (n % nrec == 0){
			printf("n = %i / %i Time step recorded. \n", n, N);
			m++;
			for (i = 0; i < J1; i++){
				H[m*J1+ i] = h1[i];
				Q[m*J1+ i] = q1[i];
				E[m*J1+ i] = e1[i];
				P[m*J1+ i] = p1[i];
				F[m*J1+ i] = f1[i];
			}
		}

		
		
	}
	
	// free memory
	free(h0);
	free(e0);
	free(f0);
	free(q0);
	free(h1);
	free(e1);
	free(f1);
	free(q1);
	free(p0);
	free(p1);

}




//// back-calculate p, q, and vs
//void fdaux(int J, int N, int M, double dr, double dt, double *params,
//           double *H, double *E, double *F, double *Q,
//		 double *gamma1, double *gamma2, double *h1, double *h2, double *pdyn, double *v1, double *v2){
//	int i, j, n, m, nrec;
//	int ii, jj;
//	int JJ = J*J;
//	int J1 = J+1;
//	int N1 = N+1;
//	int M1 = M+1;
//	int J1J1 = J1*J1;
//	double dr2 = dr*dr;
//	double dth = 0.5*dt;
//	double temp;
//	
//
//	// allocate memory and initialize to zero
//
//	// operator matrices
//	double *U      = (double*) calloc(JJ, sizeof(double));
//	double *DF     = (double*) calloc(JJ, sizeof(double));
//	double *DB     = (double*) calloc(JJ, sizeof(double));
//	double *DC     = (double*) calloc(JJ, sizeof(double));
//	double *DA     = (double*) calloc(JJ, sizeof(double));
//	double *DD     = (double*) calloc(JJ, sizeof(double));
//	double *DDDD   = (double*) calloc(JJ, sizeof(double));
//	double *DP     = (double*) calloc(JJ, sizeof(double));
//	double *DN     = (double*) calloc(JJ, sizeof(double));
//
//	double *hn     = (double*) calloc(J1, sizeof(double));
//	double *en     = (double*) calloc(J1, sizeof(double));
//	double *fn     = (double*) calloc(J1, sizeof(double));
//	double *qn     = (double*) calloc(J1, sizeof(double));
//	double *hn1    = (double*) calloc(J1, sizeof(double));
//	double *en1    = (double*) calloc(J1, sizeof(double));
//	double *fn1    = (double*) calloc(J1, sizeof(double));
//	double *qn1    = (double*) calloc(J1, sizeof(double));
//	double *hsum   = (double*) calloc(J, sizeof(double));
//	double *hdiff   = (double*) calloc(J, sizeof(double));
//	double *hdiff_rhs   = (double*) calloc(J, sizeof(double));
//
//	double *vs     = (double*) calloc(J, sizeof(double));
//	double *vsbar  = (double*) calloc(J, sizeof(double));
//	double *d1     = (double*) calloc(J, sizeof(double));
//	double *d2     = (double*) calloc(J, sizeof(double));
//	double *A      = (double*) calloc(JJ, sizeof(double));
//	double *B      = (double*) calloc(JJ, sizeof(double));
//	double *diagE  = (double*) calloc(JJ, sizeof(double));
//	double *diagF  = (double*) calloc(JJ, sizeof(double));
//
//	double *vsm    = (double*) calloc(J, sizeof(double));
//	double *hm     = (double*) calloc(J, sizeof(double));
//	double *qm     = (double*) calloc(J, sizeof(double));
//	
//	// initialize record counter and wait time
//	m = -1;
//	nrec = N/M;
//
//	// physical parameters
//	double Ca    = params[0]; double k0 = 1.0/Ca;
//	double Bo    = params[1]; double Fr = Ca/Bo;
//	double Ma    = params[2];
//	double tstop = params[3];
//	double R1    = params[4];
//	double po = params[9];
//	double pobar = params[10];
//	
//	// operator matrices
//	fdmat(J, dr, U, DF, DB, DC, DA, DD, DP, DN, DDDD);
//
//	// variable conversion
//	for (n = 0; n < N; n++){
//		if (n % nrec == 0){
//			m++;
//			// compute gamma1, gamma2, and hsum
//			for (i = 0; i < J1; i++){
//				// call solution vectors
//				hn[i]  = H[    m*J1 + i];
//				en[i]  = E[    m*J1 + i];
//				qn[i]  = Q[    m*J1 + i];
//				if (m < M) {
//					hn1[i] = H[(m+1)*J1+ i];
//					en1[i] = E[(m+1)*J1+ i];
//					qn1[i] = Q[(m+1)*J1+ i];
//				}else{
//					hn1[i] = hn[i];
//					en1[i] = en[i];
//					qn1[i] = qn[i];
//				}
//
//				// compute gamma1 and gamma2
//				gamma1[m*J1+i] = (en[i] + fn[i])/2.0;
//				gamma2[m*J1+i] = (en[i] + fn[i])/2.0;
//
//			}
//
//			// compute hsum, hdiff, h1, and h2
//			for (i = 0; i < J; i++){
//				hsum[i]	= hn[i] + (1.0 + i*i*dr*dr/2.0-dt*n);
//				hdiff_rhs[i] = dr*dr*(Bo*hsum[i]+2*pobar);
//			}
//			latri(J, DD, hdiff_rhs, hdiff);
//			for (i = 0; i < J; i++){
//				h1[m*J1+i] = (hsum[i]+hdiff[i])/2;
//				h2[m*J1+i] = (hsum[i]-hdiff[i])/2;
//			}
//			h1[m*J1+J] = 0.0;
//			h2[m*J1+J] = 1-R1*R1/2;
//
//			
//
//
//			// compute vs and vsbar then convert to v1 and v2
//			// assemble diagonal matrices
//			for (i = 0; i < J; i++){
//				diagE[i*J + i] = en[i]+en1[i]+4;
//				diagF[i*J + i] = fn[i]+fn1[i];
//				d1[i]		= 2/(nrec*dt)*(en[i]-en1[i]);
//				d2[i]		= 2/(nrec*dt)*(fn[i]-fn1[i]);
//			}
//			laband(J, DA, 0, 1, diagE, 0, 0, A);
//			laband(J, DA, 0, 1, diagF, 0, 0, B);
//			// compute vs and vsbar
//			labitri(J, A, B, B, A, d1, d2, vs, vsbar);
//			// convert vs and vsbar to v1 and v2
//			for (i = 1; i < J1; i++){
//				v1[m*J1+i] = (vs[i] + vsbar[i])*0.5;
//				v2[m*J1+i] = (vs[i] - vsbar[i])*0.5;
//			}
//			v1[m*J1+0] = 0.0;
//			v2[m*J1+0] = 0.0;
//
//			// compute dynamic pressure
//			for(i = 0; i < J; i++){
//				hm[i] = (hn[i]+hn[i+1])*0.5;
//				qm[i] = (qn[i]+qn[i+1])*0.5;
//				vsm[i] = (v1[i]+v2[i]+v1[i+1]+v2[i+1])*0.5;
//			}
//
//
//
//
//			pdyn[ m*J1 + J ] = 0.0;
//			for (i = J-1; i > -1; i--) {
//				temp = (qm[i]-hm[i]*vsm[i]/2)/hm[i]/hm[i]/hm[i];
//				pdyn[ m*J + i ] = pdyn[ m*J + i+1 ] - 12*Ca*dr*temp;
//			}
//
//			
//		}
//	
//	}
//	
//	
//	// free memory
//	free(U    );
//	free(DP   );
//	free(DN   );
//	free(DF   );
//	free(DB   );
//	free(DC   );
//	free(DA   );
//	free(DD   );
//	free(DDDD );
//	free(hn   );
//	free(en   );
//	free(fn   );
//	free(qn   );
//	free(hn1  );
//	free(en1  );
//	free(fn1  );
//	free(qn1  );
//	free(vs   );
//	free(vsbar);
//	free(d1   );
//	free(d2   );
//	free(A    );
//	free(B    );
//	free(diagE);
//	free(diagF);
//	free(hm);
//	free(qm);
//	free(vsm);
//	free(hsum);
//	free(hdiff);
//	free(hdiff_rhs);
//
//
//}

#endif

