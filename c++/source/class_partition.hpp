//----------------------------------------------------------------------------//
/*
  File:         class_partition.hpp

  Description:  Header file for class_partition
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Preprocessor def
//----------------------------------------------------------------------------//
#ifndef __incl_class_partition
#define __incl_class_partition

//----------------------------------------------------------------------------//
// Include header files
//----------------------------------------------------------------------------//
// System
#include <stdlib.h>
#include <iostream>
#include <vector>
// External
#include "mpi.h"
// Internal

//----------------------------------------------------------------------------//
// Partition information class
//----------------------------------------------------------------------------//
class class_partition
{
	private:
//----------------------------------------------------------------------------//
// Private variables
//----------------------------------------------------------------------------//

	public:
//----------------------------------------------------------------------------//
// Public variables
//----------------------------------------------------------------------------//
	int ncell[3];
	int icell[3];
	double dx[3];
};

//----------------------------------------------------------------------------//
// Prototype module subroutines
//----------------------------------------------------------------------------//
	std::vector<class_partition> partition_setup( int, int [3], int [3], double [3], bool, bool, bool );

#endif
