//----------------------------------------------------------------------------//
/*
  File:         setup2d.cpp

  Description:  Sets up the greenfish object for 2D simulations
*/
//----------------------------------------------------------------------------//

void class_greenfish::setup2d( int ncell[2], int bound_cond[2], double dx[2] )
{

//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int nproc, rank;

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//----------------------------------------------------------------------------//
// Pass domain info to partitions
//----------------------------------------------------------------------------//
	domain_ncell[0] = ncell[0];
	domain_ncell[1] = ncell[1];
	domain_ncell[2] = 1;

	domain_bc[0] = bound_cond[0];
	domain_bc[1] = bound_cond[1];
	domain_bc[2] = -1;

	domain_dx[0] = dx[0];
	domain_dx[1] = dx[1];
	domain_dx[2] = 0.0;

//----------------------------------------------------------------------------//
// Setup real topology
//----------------------------------------------------------------------------//
	partition = partition_setup( 0, domain_ncell, domain_bc, domain_dx, false, false, false );

//----------------------------------------------------------------------------//
// Setup Greens function
//----------------------------------------------------------------------------//
	greens2d( domain_ncell, domain_bc, domain_dx, regularisation );

//----------------------------------------------------------------------------//
// Setup pencil topologies
//----------------------------------------------------------------------------//
	xpen = partition_setup( 0, domain_ncell, domain_bc, domain_dx, true, false, false );
	ypen = partition_setup( 1, domain_ncell, domain_bc, domain_dx, true, false, false );

//----------------------------------------------------------------------------//
// Construct communicators between topologies
//----------------------------------------------------------------------------//
	real2xpen = class_communication( partition, xpen );
	xpen2real = class_communication( xpen, partition );

	xpen2ypen = class_communication( xpen, ypen );
	ypen2xpen = class_communication( ypen, xpen );

	return;
}

