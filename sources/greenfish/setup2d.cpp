//----------------------------------------------------------------------------//
/*
  File:         setup2d.cpp

  Description:  
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
	partition.real = partition_setup2d( 0, ncell, bound_cond, dx, false, false );

//----------------------------------------------------------------------------//
// Setup Greens function
//----------------------------------------------------------------------------//
	greens2d();

//----------------------------------------------------------------------------//
// Setup pencil topologies
//----------------------------------------------------------------------------//
	partition.xpen = partition_setup2d( 0, ncell, bound_cond, dx, true, false );
	partition.ypen = partition_setup2d( 1, ncell, bound_cond, dx, true, false );

//----------------------------------------------------------------------------//
// Construct communicators between topologies
//----------------------------------------------------------------------------//
	real2xpen.setup2d( partition.real, partition.xpen );
	xpen2real.setup2d( partition.xpen, partition.real );

	xpen2ypen.setup2d( partition.xpen, partition.ypen );
	ypen2xpen.setup2d( partition.ypen, partition.xpen );

	return;
}

