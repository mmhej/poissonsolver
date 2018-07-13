//----------------------------------------------------------------------------//
/*
  File:         class_greenfpm.hpp

  Description:  Header file for class_greenfpm.cpp
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Preprocessor def
//----------------------------------------------------------------------------//
#ifndef __incl_class_greenfish
#define __incl_class_greenfish

//----------------------------------------------------------------------------//
// Include header files
//----------------------------------------------------------------------------//
// System
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
class class_greenfish
{

	private:
//----------------------------------------------------------------------------//
// Private variables
//----------------------------------------------------------------------------//
		int    domain_ncell[3];
		int    domain_bc[3];
		double domain_dx[3];

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

		std::complex<double> * mapG = NULL;
		std::complex<double> * rhsG = NULL;

//----------------------------------------------------------------------------//
// Prototype private routines
//----------------------------------------------------------------------------//
		void map2d( class_communication );
		void greens2d( );

	public:
//----------------------------------------------------------------------------//
// Public variables
//----------------------------------------------------------------------------//
		class_partition partition;

//----------------------------------------------------------------------------//
// Prototype public routines
//----------------------------------------------------------------------------//
		void setup2d( int [2], int [2], double [2] );
		void solve2d(  );

		void push2d( double *, double *, double *, double *, double *, double * );
		void pull2d( double *, double *, double *, double *, double *, double * );

};


#endif
