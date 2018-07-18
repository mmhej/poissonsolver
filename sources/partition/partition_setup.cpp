//----------------------------------------------------------------------------//
/*
  File:         partition_setup.cpp

  Description:  
*/
//----------------------------------------------------------------------------//

std::vector<class_partition> partition_setup( int pencil_dir, 
                                              int domain_ncell[3], 
                                              int bound_cond[3],
                                              double dx[3],
                                              bool extendX,
                                              bool extendY,
                                              bool extendZ )
{

//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int i,j,k;
	int nproc, rank;
	int nsub[3];
	int dom_ncell[3];
	int ncell[3];
	int icell[3];

	int rem_ncell[3];
	int iproc;

//----------------------------------------------------------------------------//
// Objects
//----------------------------------------------------------------------------//
	class_partition part;
	std::vector<class_partition> part_all;

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
		nsub[2] = 1;
	}
	else if(pencil_dir == 1)
	{
		nsub[0] = nproc;
		nsub[1] = 1;
		nsub[2] = 1;
	}
	else
	{
		std::cerr << " [greenfish.partition_setup3d]: Error." << std::endl;
		exit(EXIT_FAILURE);
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

	if(bound_cond[2] == 0 && extendZ)
	{
		dom_ncell[2] = 2*domain_ncell[2];
	}
	else
	{
		dom_ncell[2] = domain_ncell[2];
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

			icell[2] = 0;
			for(k = 0; k < nsub[2]; ++k)
			{


//----------------------------------------------------------------------------//
// Divide the grid points to subdomain (ghost not included)
// If the grid points does not equally divide with the number of processors
// the remaining grid points are distributed equally over the last processors
//----------------------------------------------------------------------------//
				rem_ncell[0] = dom_ncell[0] % nsub[0];
				rem_ncell[1] = dom_ncell[1] % nsub[1];
				rem_ncell[2] = dom_ncell[2] % nsub[2];

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

				if(k >= nsub[2]-rem_ncell[2])
				{
					ncell[2] = dom_ncell[2]/nsub[2] + 1;
				}
				else
				{
					ncell[2] = dom_ncell[2]/nsub[2];
				}

//----------------------------------------------------------------------------//
// Store the partition
//----------------------------------------------------------------------------//
				part.ncell[0] = ncell[0];
				part.ncell[1] = ncell[1];
				part.ncell[2] = ncell[2];

				part.icell[0] = icell[0];
				part.icell[1] = icell[1];
				part.icell[2] = icell[2];

				part.dx[0]    = dx[0];
				part.dx[1]    = dx[1];
				part.dx[2]    = dx[2];

				part_all.push_back( part );

//----------------------------------------------------------------------------//
// Advance processor and cell counters
//----------------------------------------------------------------------------//
				++iproc;
				icell[2] = icell[2] + ncell[2];
			}
			icell[1] = icell[1] + ncell[1];
		}
		icell[0] = icell[0] + ncell[0];
	}


	return part_all;
}
