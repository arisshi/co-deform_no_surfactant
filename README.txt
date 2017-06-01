/*-------------------------------------------------------------------------------*/
/*--     THIN-FILM EVOLUTION OVER DEFORMABLE SPHERE W/O SURFACTANT (V13)       --*/
/*-------------------------------------------------------------------------------*/

Description: 
  	Numerical solution of the coupled, one-dimensional evolution of film height over a deformable bubble. 
  	Algorithm implements the finite difference method / quasi-linearization to solve a system
	of one first order hyperbolic PDE and two second-order parabolic PDEs.
	See supporting documentation for the derivation of the basic equations and
	solution algorithm. CODE GIVES NON-PHYSICAL RESULTS.

Language: C++

Libraries: LAPACK, ATLAS, BLAS, GSL, EIGEN

Directory tree:
	include - contains header files
	src - contains input parameters and executable
	src/output - stores output files written by simulation (in .txt format)
	postproc - contains scripts for post-processing	
	
Implementations (header files):
	fd.h - finite difference subroutines
	la.h - linear algebra subroutines
	rw.h - read / write subroutines
	bessel.h - Bessel function evaluation subroutine (written by Kapteyn Institute Groningen)

Input parameters (in src/params.in)
	CA - capillary number
	BO - bond number
	MA - Marangoni number (Set to zero for this no surfactant system)
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

Authors:
	Joseph Barakat and Xingyi Shi
