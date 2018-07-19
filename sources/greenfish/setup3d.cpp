//----------------------------------------------------------------------------//
/*
  File:         setup3d.cpp

  Description:  
*/
//----------------------------------------------------------------------------//

void class_greenfish::setup3d( int ncell[3], int bound_cond[3], double dx[3] )
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
	domain_ncell[2] = ncell[2];

	domain_bc[0] = bound_cond[0];
	domain_bc[1] = bound_cond[1];
	domain_bc[2] = bound_cond[2];

	domain_dx[0] = dx[0];
	domain_dx[1] = dx[1];
	domain_dx[2] = dx[2];

//----------------------------------------------------------------------------//
// Setup real topology
//----------------------------------------------------------------------------//
	partition = partition_setup( 0, ncell, bound_cond, dx, false, false, false );

//----------------------------------------------------------------------------//
// Setup Greens function
//----------------------------------------------------------------------------//
	greens3d();

//----------------------------------------------------------------------------//
// Setup pencil topologies
//----------------------------------------------------------------------------//
	xpen = partition_setup( 0, domain_ncell, domain_bc, domain_dx, true, false, false );
	ypen = partition_setup( 1, domain_ncell, domain_bc, domain_dx, true,  true, false );

//----------------------------------------------------------------------------//
// Construct communicators between topologies
//----------------------------------------------------------------------------//
	real2xpen = class_communication( partition, xpen );
	xpen2real = class_communication( xpen, partition );

	xpen2ypen = class_communication( xpen, ypen );
	ypen2xpen = class_communication( ypen, xpen );

	return;
}

