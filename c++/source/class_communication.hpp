//----------------------------------------------------------------------------//
/*
  File:         class_communication.hpp

  Description:  Header file for class_communication
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Preprocessor def
//----------------------------------------------------------------------------//
#ifndef __incl_class_communication
#define __incl_class_communication

//----------------------------------------------------------------------------//
// Include header files
//----------------------------------------------------------------------------//
// System
#include <iostream>
#include <vector>
// External
#include "mpi.h"
// Internal
#include "class_partition.hpp"


//----------------------------------------------------------------------------//
// Communication information class
//----------------------------------------------------------------------------//
class class_cominfo
{
	private:
//----------------------------------------------------------------------------//
// Private variables
//----------------------------------------------------------------------------//

	public:
//----------------------------------------------------------------------------//
// Public variables
//----------------------------------------------------------------------------//
		int iproc;
		int jproc;

		int nway;

		int i2j_ncell[3];
		int i2j_min_send[3];
		int i2j_max_send[3];
		int i2j_min_recv[3];
		int i2j_max_recv[3];

		int i2j_ncell_partition_send[3];
		int i2j_ncell_partition_recv[3];

		int j2i_ncell[3];
		int j2i_min_send[3];
		int j2i_max_send[3];
		int j2i_min_recv[3];
		int j2i_max_recv[3];

		int j2i_ncell_partition_send[3];
		int j2i_ncell_partition_recv[3];

};


class class_partinfo
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

};



//----------------------------------------------------------------------------//
// Communication class
//----------------------------------------------------------------------------//
class class_communication
{
	private:
//----------------------------------------------------------------------------//
// Private variables
//----------------------------------------------------------------------------//


	public:
//----------------------------------------------------------------------------//
// Public variables
//----------------------------------------------------------------------------//
		int ncomm;
		std::vector< class_partinfo > partition_send;
		std::vector< class_partinfo > partition_recv;

		std::vector< class_cominfo > info;

//----------------------------------------------------------------------------//
// Prototype public subroutines
//----------------------------------------------------------------------------//
		class_communication( void );
		class_communication( std::vector<class_partition>, 
		                     std::vector<class_partition> );


};


#endif
