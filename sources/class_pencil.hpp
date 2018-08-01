//----------------------------------------------------------------------------//
/*
  File:         class_pencil.hpp

  Description:  Header file for class_pencil
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Preprocessor def
//----------------------------------------------------------------------------//
#ifndef __incl_class_pencil
#define __incl_class_pencil

//----------------------------------------------------------------------------//
// Include header files
//----------------------------------------------------------------------------//
// System
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
// External
#include "mpi.h"
// Internal
//#include "class_topology.hpp"


//----------------------------------------------------------------------------//
// Pencil class
//----------------------------------------------------------------------------//
class class_pencil
{
	private:
//----------------------------------------------------------------------------//
// Private variables
//----------------------------------------------------------------------------//

	public:
//----------------------------------------------------------------------------//
// Public variables
//----------------------------------------------------------------------------//
		int nfft;

		bool bX;
		bool bY;
		bool bZ;

		std::complex<double> * X    = NULL;
		std::complex<double> * Y    = NULL;
		std::complex<double> * Z    = NULL;

		std::complex<double> * auxX = NULL;
		std::complex<double> * auxY = NULL;
		std::complex<double> * auxZ = NULL;

//----------------------------------------------------------------------------//
// Prototype subroutines
//----------------------------------------------------------------------------//
		class_pencil( bool, bool, bool );
		~class_pencil( void );

		void resize( int );

		void fft_shift( void );

		void fft( void );
		void ifft( void );

};


#endif
