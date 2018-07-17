//----------------------------------------------------------------------------//
/*
  File:         partition_setup2d.cpp

  Description:  
*/
//----------------------------------------------------------------------------//

std::vector<class_partition_info> partition_setup2d( int pencil_dir, 
                                                     int domain_ncell[2], 
                                                     int bound_cond[2],
                                                     double dx[2],
                                                     bool extendX,
                                                     bool extendY )
{

//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int i,j;
	int nproc, rank;
	int nsub[2];
	int dom_ncell[2];
	int ncell[2];
	int icell[2];

	int rem_ncell[2];
	int iproc;

//----------------------------------------------------------------------------//
// Objects
//----------------------------------------------------------------------------//
	class_partition_info partition;
	std::vector<class_partition_info> partition_all;

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//============================================================================//
// Setup pencil partition
//============================================================================//
//----------------------------------------------------------------------------//
// Partition into pencil partition given the input direction
//----------------------------------------------------------------------------//
	if(pencil_dir == 0)
	{
		nsub[0] = 1;
		nsub[1] = nproc;
	}
	else if(pencil_dir == 1)
	{
		nsub[0] = nproc;
		nsub[1] = 1;
	}

//----------------------------------------------------------------------------//
// Extend domain in the X-direction if this is unbounded BCs (y-direction 
// extension will be handled on pencils)
//----------------------------------------------------------------------------//
	if(bound_cond[0] == 0 && extendX)
	{
		dom_ncell[0] = 2*domain_ncell[0];
	}
	else
	{
		dom_ncell[0] = domain_ncell[0];
	}

	if(bound_cond[1] == 0 && extendY)
	{
		dom_ncell[1] = 2*domain_ncell[1];
	}
	else
	{
		dom_ncell[1] = domain_ncell[1];
	}


//----------------------------------------------------------------------------//
// Assign subdomains and keep track of the local index of the subdomain on
// the processor
//----------------------------------------------------------------------------//
	iproc = 0;

	icell[0] = 0;
	for(i = 0; i < nsub[0]; ++i)
	{

		icell[1] = 0;
		for(j = 0; j < nsub[1]; ++j)
		{

//----------------------------------------------------------------------------//
// Divide the grid points to subdomain (ghost not included)
// If the grid points does not equally divide with the number of processors
// the remaining grid points are distributed equally over the last processors
//----------------------------------------------------------------------------//
			rem_ncell[0] = dom_ncell[0] % nsub[0];
			rem_ncell[1] = dom_ncell[1] % nsub[1];

			if(i >= nsub[0]-rem_ncell[0])
			{
				ncell[0] = dom_ncell[0]/nsub[0] + 1;
			}
			else
			{
				ncell[0] = dom_ncell[0]/nsub[0];
			}

			if(j >= nsub[1]-rem_ncell[1])
			{
				ncell[1] = dom_ncell[1]/nsub[1] + 1;
			}
			else
			{
				ncell[1] = dom_ncell[1]/nsub[1];
			}

//----------------------------------------------------------------------------//
// Store the partition
//----------------------------------------------------------------------------//
			partition.ncell[0] = ncell[0];
			partition.ncell[1] = ncell[1];
			partition.ncell[2] = 1;

			partition.icell[0] = icell[0];
			partition.icell[1] = icell[1];
			partition.icell[2] = 0;

			partition.dx[0]    = dx[0];
			partition.dx[1]    = dx[1];
			partition.dx[2]    = 0.0;

			partition_all.push_back( partition );

//----------------------------------------------------------------------------//
// Advance processor and cell counters
//----------------------------------------------------------------------------//
			++iproc;
			icell[1] = icell[1] + ncell[1];
		}

		icell[0] = icell[0] + ncell[0];
	}


	return partition_all;
}
