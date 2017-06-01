/*-------------------------------------------------------------------------------*/
/*--  THIN-FILM EVOLUTION OVER RIGID SPHERE W/INSOLUBLE SURFACTANT (DRAFT #06) --*/
/*-------------------------------------------------------------------------------*/

Description: 
  Numerical solution of the coupled, one-dimensional evolution of film height and 
	insoluble surfactant concentration over a rigid, spherical substrate. Algorithm 
	implements the finite difference method / quasi-linearization to solve a system
	of two fourth-order, parabolic PDEs (bipentadiagonal linear system when quasi-
	linearized).
	See supporting documentation for the derivation of the basic equations and
	solution algorithm.

Language: C++

Libraries: LAPACK, ATLAS, BLAS, GSL

Directory tree:
	include - contains header files
	src - contains input parameters and executable
	output - stores output files written by simulation (in .txt format)
	postproc - contains scripts for post-processing	
	
Implementations (header files):
	fd.h - finite difference subroutines
	la.h - linear algebra subroutines
	rw.h - read / write subroutines

Input parameters (in src/params.in)
	CA - capillary number
	BO - bond number
	MA - Marangoni number
	TSTOP - time when sphere stops moving
	R1 - upper limit of radial domain
	T1 - upper limit of time domain 
	DR - spatial grid size
	DT - time step size
	DTREC - time points to output data file

Instructions:
	1. Make executable,
		cd src
		make clean
		make
	2. Set input parameters in src/params.in
	3. Run executable,
		cd src
		./run
	4. Post-process output files written to output.

Notes:
	Since the problem is fourth-order in space, the CFL number is given by
	
	        dt D
	  CFL = ----
	        dr^4
	
	where D is a fourth-order diffusion coefficient. Typically,
	
	  D ~ O(r1^6 / Ca)    or    D ~ O(r1^4 Ma)
	
	Having CFL = O(1) is desirable for numerical accuracy, but is not 
	required for numerical stability.