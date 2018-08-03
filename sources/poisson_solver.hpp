//----------------------------------------------------------------------------//
/*
  File:         poisson_solver.hpp

  Description:  Header file for poisson_solver
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Preprocessor def
//----------------------------------------------------------------------------//
#ifndef __incl_poisson_solver
#define __incl_poisson_solver

//----------------------------------------------------------------------------//
// Include header files
//----------------------------------------------------------------------------//
// System
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
// External
#include "mpi.h"
// Internal
#include "class_partition.hpp"
#include "class_communication.hpp"
#include "class_pencil.hpp"

// outputting to file
#include <iomanip>
#include <fstream>
#include <string>

//----------------------------------------------------------------------------//
// Green library class
//----------------------------------------------------------------------------//
class poisson_solver
{

	private:
//----------------------------------------------------------------------------//
// Private variables
//----------------------------------------------------------------------------//
		int    domain_ncell[3];
		int    domain_bc[3];
		double domain_dx[3];

		std::vector<class_partition> xpen;
		std::vector<class_partition> ypen;
		std::vector<class_partition> zpen;

		class_communication real2xpen;
		class_communication xpen2real;

		class_communication xpen2ypen;
		class_communication ypen2xpen;

		class_communication xpen2zpen;
		class_communication zpen2xpen;

		class_communication ypen2zpen;
		class_communication zpen2ypen;

		std::complex<double> * rhsX = NULL;
		std::complex<double> * rhsY = NULL;
		std::complex<double> * rhsZ = NULL;

		std::complex<double> * lhsX = NULL;
		std::complex<double> * lhsY = NULL;
		std::complex<double> * lhsZ = NULL;

		std::complex<double> * ikX  = NULL;
		std::complex<double> * ikY  = NULL;
		std::complex<double> * ikZ  = NULL;

		double * zeta = NULL;

		std::complex<double> * mapG = NULL;
		std::complex<double> * G2D  = NULL;
		std::complex<double> * G3D  = NULL;

//----------------------------------------------------------------------------//
// Prototype private routines
//----------------------------------------------------------------------------//
		void map( class_communication );
		void greens2d( int [3], int [3], double [3], int );
		void greens3d( int [3], int [3], double [3], int );

	public:
//----------------------------------------------------------------------------//
// Public variables
//----------------------------------------------------------------------------//
		std::vector<class_partition> partition;

// Options
		bool lhs_grad      = false;
		bool lhs_div       = false;
		bool lhs_curl      = false;
		bool rhs_reproject = false;
		int  regularisation = 0;

//----------------------------------------------------------------------------//
// Prototype public routines
//----------------------------------------------------------------------------//
		poisson_solver( void );
		~poisson_solver( void );

		void setup2d( int [2], int [2], double [2] );
		void setup3d( int [3], int [3], double [3] );

		void solve2d(  );
		void solve3d(  );

		void push( double *, double *, double *, double *, double *, double * );
		void pull( double *, double *, double *, double *, double *, double * );

};


#endif
