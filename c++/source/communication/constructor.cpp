//----------------------------------------------------------------------------//
/*
  File:         constructor.cpp

  Description:  Constructor for class_communication
*/
//----------------------------------------------------------------------------//
class_communication::class_communication( void ){}

class_communication::class_communication( 
                                  std::vector<class_partition> partition_send, 
                                  std::vector<class_partition> partition_recv )
{

//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int i,j;
	int nproc, rank;
	int iproc, jproc;
	int icomm, ncomm, sum_comm;
	int imin[3], imax[3], jmin[3], jmax[3];
	int ncell_send[3], ncell_recv[3];

//----------------------------------------------------------------------------//
// Objects
//----------------------------------------------------------------------------//
	class_cominfo tmp_info;
	class_partinfo prt_info;
	std::vector< class_cominfo > tmp_info_all;

	std::vector< int > proc_busy;
	std::vector< int > comm_list;

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//----------------------------------------------------------------------------//
// Store partition info
//----------------------------------------------------------------------------//
	for(i = 0; i < nproc; ++i)
	{
		prt_info.ncell[0] = partition_send[i].ncell[0];
		prt_info.ncell[1] = partition_send[i].ncell[1];
		prt_info.ncell[2] = partition_send[i].ncell[2];

		class_communication::partition_send.push_back( prt_info );

		prt_info.ncell[0] = partition_recv[i].ncell[0];
		prt_info.ncell[1] = partition_recv[i].ncell[1];
		prt_info.ncell[2] = partition_recv[i].ncell[2];

		class_communication::partition_recv.push_back( prt_info );
	}

//----------------------------------------------------------------------------//
// Find number of communications
//----------------------------------------------------------------------------//
	icomm = 0;
	for( iproc = 0; iproc < nproc; ++iproc )
	{
		for( jproc = iproc; jproc < nproc; ++jproc )
		{

//----------------------------------------------------------------------------//
// Processors and number of communication ways
//----------------------------------------------------------------------------//
			tmp_info.iproc = iproc;
			tmp_info.jproc = jproc;
			tmp_info.nway = 0;

//----------------------------------------------------------------------------//
// Overlap of subdomains from iproc to jproc
//----------------------------------------------------------------------------//
			imin[0] = partition_send[iproc].icell[0];
			imin[1] = partition_send[iproc].icell[1];
			imin[2] = partition_send[iproc].icell[2];

			imax[0] = partition_send[iproc].icell[0] + partition_send[iproc].ncell[0] - 1;
			imax[1] = partition_send[iproc].icell[1] + partition_send[iproc].ncell[1] - 1;
			imax[2] = partition_send[iproc].icell[2] + partition_send[iproc].ncell[2] - 1;

			jmin[0] = partition_recv[jproc].icell[0];
			jmin[1] = partition_recv[jproc].icell[1];
			jmin[2] = partition_recv[jproc].icell[2];

			jmax[0] = partition_recv[jproc].icell[0] + partition_recv[jproc].ncell[0] - 1;
			jmax[1] = partition_recv[jproc].icell[1] + partition_recv[jproc].ncell[1] - 1;
			jmax[2] = partition_recv[jproc].icell[2] + partition_recv[jproc].ncell[2] - 1;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Check overlap
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			if( jmax[0] >= imin[0] && jmax[1] >= imin[1] && jmax[2] >= imin[2] &&
			    imax[0] >= jmin[0] && imax[1] >= jmin[1] && imax[2] >= jmin[2] )
			{

				tmp_info.nway += 1;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Overlaping indexes - min index
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				if( imin[0] == jmin[0] )
				{
					tmp_info.i2j_min_send[0] = 0;
					tmp_info.i2j_min_recv[0] = 0;
				}
				else if( imin[0] > jmin[0] )
				{
					tmp_info.i2j_min_send[0] = 0;
					tmp_info.i2j_min_recv[0] = imin[0] - jmin[0];
				}
				else if( imin[0] < jmin[0] )
				{
					tmp_info.i2j_min_send[0] = jmin[0] - imin[0];
					tmp_info.i2j_min_recv[0] = 0;
				}

				if( imin[1] == jmin[1] )
				{
					tmp_info.i2j_min_send[1] = 0;
					tmp_info.i2j_min_recv[1] = 0;
				}
				else if( imin[1] > jmin[1] )
				{
					tmp_info.i2j_min_send[1] = 0;
					tmp_info.i2j_min_recv[1] = imin[1] - jmin[1];
				}
				else if( imin[1] < jmin[1] )
				{
					tmp_info.i2j_min_send[1] = jmin[1] - imin[1];
					tmp_info.i2j_min_recv[1] = 0;
				}

				if( imin[2] == jmin[2] )
				{
					tmp_info.i2j_min_send[2] = 0;
					tmp_info.i2j_min_recv[2] = 0;
				}
				else if( imin[2] > jmin[2] )
				{
					tmp_info.i2j_min_send[2] = 0;
					tmp_info.i2j_min_recv[2] = imin[2] - jmin[2];
				}
				else if( imin[2] < jmin[2] )
				{
					tmp_info.i2j_min_send[2] = jmin[2] - imin[2];
					tmp_info.i2j_min_recv[2] = 0;
				}


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Overlaping indexes - max index
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				if( imax[0] == jmax[0] )
				{
					tmp_info.i2j_max_send[0] = partition_send[iproc].ncell[0] - 1;
					tmp_info.i2j_max_recv[0] = partition_recv[jproc].ncell[0] - 1;
				}
				else if( imax[0] < jmax[0] )
				{
					tmp_info.i2j_max_send[0] = partition_send[iproc].ncell[0] - 1;
					tmp_info.i2j_max_recv[0] = partition_recv[jproc].ncell[0]
					                         - 1 - (jmax[0] - imax[0]);
				}
				else if( imax[0] > jmax[0] )
				{
					tmp_info.i2j_max_send[0] = partition_send[iproc].ncell[0]
					                         - 1 - (imax[0] - jmax[0]);
					tmp_info.i2j_max_recv[0] = partition_recv[jproc].ncell[0] - 1;
				}

				if( imax[1] == jmax[1] )
				{
					tmp_info.i2j_max_send[1] = partition_send[iproc].ncell[1] - 1;
					tmp_info.i2j_max_recv[1] = partition_recv[jproc].ncell[1] - 1;
				}
				else if( imax[1] < jmax[1] )
				{
					tmp_info.i2j_max_send[1] = partition_send[iproc].ncell[1] - 1;
					tmp_info.i2j_max_recv[1] = partition_recv[jproc].ncell[1]
					                         - 1 - (jmax[1] - imax[1]);
				}
				else if( imax[1] > jmax[1] )
				{
					tmp_info.i2j_max_send[1] = partition_send[iproc].ncell[1]
					                         - 1 - (imax[1] - jmax[1]);
					tmp_info.i2j_max_recv[1] = partition_recv[jproc].ncell[1] - 1;
				}

				if( imax[2] == jmax[2] )
				{
					tmp_info.i2j_max_send[2] = partition_send[iproc].ncell[2] - 1;
					tmp_info.i2j_max_recv[2] = partition_recv[jproc].ncell[2] - 1;
				}
				else if( imax[2] < jmax[2] )
				{
					tmp_info.i2j_max_send[2] = partition_send[iproc].ncell[2] - 1;
					tmp_info.i2j_max_recv[2] = partition_recv[jproc].ncell[2]
					                         - 1 - (jmax[2] - imax[2]);
				}
				else if( imax[2] > jmax[2] )
				{
					tmp_info.i2j_max_send[2] = partition_send[iproc].ncell[2]
					                         - 1 - (imax[2] - jmax[2]);
					tmp_info.i2j_max_recv[2] = partition_recv[jproc].ncell[2] - 1;
				}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store number of cells to be communicated
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				ncell_send[0] = tmp_info.i2j_max_send[0] - tmp_info.i2j_min_send[0] + 1;
				ncell_send[1] = tmp_info.i2j_max_send[1] - tmp_info.i2j_min_send[1] + 1;
				ncell_send[2] = tmp_info.i2j_max_send[2] - tmp_info.i2j_min_send[2] + 1;

				ncell_recv[0] = tmp_info.i2j_max_recv[0] - tmp_info.i2j_min_recv[0] + 1;
				ncell_recv[1] = tmp_info.i2j_max_recv[1] - tmp_info.i2j_min_recv[1] + 1;
				ncell_recv[2] = tmp_info.i2j_max_recv[2] - tmp_info.i2j_min_recv[2] + 1;

				if( ncell_send[0] != ncell_recv[0] 
				 || ncell_send[1] != ncell_recv[1]
				 || ncell_send[2] != ncell_recv[2] )
				{
					std::cout << "ERROR!! number of send and recieve cells does not match" << std::endl;
				}

				tmp_info.i2j_ncell[0] = ncell_send[0];
				tmp_info.i2j_ncell[1] = ncell_send[1];
				tmp_info.i2j_ncell[2] = ncell_send[2];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store total number of cells of the two partition
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				tmp_info.i2j_ncell_partition_send[0] = partition_send[iproc].ncell[0];
				tmp_info.i2j_ncell_partition_send[1] = partition_send[iproc].ncell[1];
				tmp_info.i2j_ncell_partition_send[2] = partition_send[iproc].ncell[2];

				tmp_info.i2j_ncell_partition_recv[0] = partition_recv[jproc].ncell[0];
				tmp_info.i2j_ncell_partition_recv[1] = partition_recv[jproc].ncell[1];
				tmp_info.i2j_ncell_partition_recv[2] = partition_recv[jproc].ncell[2];

			}

//----------------------------------------------------------------------------//
// Overlap of subdomains from jproc to iproc
//----------------------------------------------------------------------------//
			jmin[0] = partition_send[jproc].icell[0];
			jmin[1] = partition_send[jproc].icell[1];
			jmin[2] = partition_send[jproc].icell[2];

			jmax[0] = partition_send[jproc].icell[0] 
			        + partition_send[jproc].ncell[0] - 1;
			jmax[1] = partition_send[jproc].icell[1] 
			        + partition_send[jproc].ncell[1] - 1;
			jmax[2] = partition_send[jproc].icell[2] 
			        + partition_send[jproc].ncell[2] - 1;

			imin[0] = partition_recv[iproc].icell[0];
			imin[1] = partition_recv[iproc].icell[1];
			imin[2] = partition_recv[iproc].icell[2];

			imax[0] = partition_recv[iproc].icell[0] + partition_recv[iproc].ncell[0] - 1;
			imax[1] = partition_recv[iproc].icell[1] + partition_recv[iproc].ncell[1] - 1;
			imax[2] = partition_recv[iproc].icell[2] + partition_recv[iproc].ncell[2] - 1;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Check overlap
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			if( imax[0] >= jmin[0] && imax[1] >= jmin[1] && imax[2] >= jmin[2] &&
			    jmax[0] >= imin[0] && jmax[1] >= imin[1] && jmax[2] >= imin[2] &&
			    iproc != jproc )
			{

				tmp_info.nway += 2;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Overlaping indexes - min index
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				if( jmin[0] == imin[0] )
				{
					tmp_info.j2i_min_send[0] = 0;
					tmp_info.j2i_min_recv[0] = 0;
				}
				else if( jmin[0] > imin[0] )
				{
					tmp_info.j2i_min_send[0] = 0;
					tmp_info.j2i_min_recv[0] = jmin[0] - imin[0];
				}
				else if( jmin[0] < imin[0] )
				{
					tmp_info.j2i_min_send[0] = imin[0] - jmin[0];
					tmp_info.j2i_min_recv[0] = 0;
				}

				if( jmin[1] == imin[1] )
				{
					tmp_info.j2i_min_send[1] = 0;
					tmp_info.j2i_min_recv[1] = 0;
				}
				else if( jmin[1] > imin[1] )
				{
					tmp_info.j2i_min_send[1] = 0;
					tmp_info.j2i_min_recv[1] = jmin[1] - imin[1];
				}
				else if( jmin[1] < imin[1] )
				{
					tmp_info.j2i_min_send[1] = imin[1] - jmin[1];
					tmp_info.j2i_min_recv[1] = 0;
				}

				if( jmin[2] == imin[2] )
				{
					tmp_info.j2i_min_send[2] = 0;
					tmp_info.j2i_min_recv[2] = 0;
				}
				else if( jmin[2] > imin[2] )
				{
					tmp_info.j2i_min_send[2] = 0;
					tmp_info.j2i_min_recv[2] = jmin[2] - imin[2];
				}
				else if( jmin[2] < imin[2] )
				{
					tmp_info.j2i_min_send[2] = imin[2] - jmin[2];
					tmp_info.j2i_min_recv[2] = 0;
				}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Overlaping indexes - max index
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				if( jmax[0] == imax[0] )
				{
					tmp_info.j2i_max_send[0] = partition_send[jproc].ncell[0] - 1;
					tmp_info.j2i_max_recv[0] = partition_recv[iproc].ncell[0] - 1;
				}
				else if( jmax[0] < imax[0] )
				{
					tmp_info.j2i_max_send[0] = partition_send[jproc].ncell[0] - 1;
					tmp_info.j2i_max_recv[0] = partition_recv[iproc].ncell[0]
					                         - 1 - (imax[0] - jmax[0]);
				}
				else if( jmax[0] > imax[0] )
				{
					tmp_info.j2i_max_send[0] = partition_send[jproc].ncell[0]
					                         - 1 - (jmax[0] - imax[0]);
					tmp_info.j2i_max_recv[0] = partition_recv[iproc].ncell[0] - 1;
				}

				if( jmax[1] == imax[1] )
				{
					tmp_info.j2i_max_send[1] = partition_send[jproc].ncell[1] - 1;
					tmp_info.j2i_max_recv[1] = partition_recv[iproc].ncell[1] - 1;
				}
				else if( jmax[1] < imax[1] )
				{
					tmp_info.j2i_max_send[1] = partition_send[jproc].ncell[1] - 1;
					tmp_info.j2i_max_recv[1] = partition_recv[iproc].ncell[1]
					                         - 1 - (imax[1] - jmax[1]);
				}
				else if( jmax[1] > imax[1] )
				{
					tmp_info.j2i_max_send[1] = partition_send[jproc].ncell[1]
					                         - 1 - (jmax[1] - imax[1]);
					tmp_info.j2i_max_recv[1] = partition_recv[iproc].ncell[1] - 1;
				}

				if( jmax[2] == imax[2] )
				{
					tmp_info.j2i_max_send[2] = partition_send[jproc].ncell[2] - 1;
					tmp_info.j2i_max_recv[2] = partition_recv[iproc].ncell[2] - 1;
				}
				else if( jmax[2] < imax[2] )
				{
					tmp_info.j2i_max_send[2] = partition_send[jproc].ncell[2] - 1;
					tmp_info.j2i_max_recv[2] = partition_recv[iproc].ncell[2]
					                         - 1 - (imax[2] - jmax[2]);
				}
				else if( jmax[2] > imax[2] )
				{
					tmp_info.j2i_max_send[2] = partition_send[jproc].ncell[2]
					                         - 1 - (jmax[2] - imax[2]);
					tmp_info.j2i_max_recv[2] = partition_recv[iproc].ncell[2] - 1;
				}


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store number of cells to be communicated
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				ncell_send[0] = tmp_info.j2i_max_send[0] - tmp_info.j2i_min_send[0] + 1;
				ncell_send[1] = tmp_info.j2i_max_send[1] - tmp_info.j2i_min_send[1] + 1;
				ncell_send[2] = tmp_info.j2i_max_send[2] - tmp_info.j2i_min_send[2] + 1;

				ncell_recv[0] = tmp_info.j2i_max_recv[0] - tmp_info.j2i_min_recv[0] + 1;
				ncell_recv[1] = tmp_info.j2i_max_recv[1] - tmp_info.j2i_min_recv[1] + 1;
				ncell_recv[2] = tmp_info.j2i_max_recv[2] - tmp_info.j2i_min_recv[2] + 1;

				if( ncell_send[0] != ncell_recv[0]
				 || ncell_send[1] != ncell_recv[1]
				 || ncell_send[2] != ncell_recv[2] )
				{
					std::cout << "ERROR!! number of send and recieve cells does not match" << std::endl;
				}

				tmp_info.j2i_ncell[0] = ncell_send[0];
				tmp_info.j2i_ncell[1] = ncell_send[1];
				tmp_info.j2i_ncell[2] = ncell_send[2];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store total number of cells of the two partition
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				tmp_info.j2i_ncell_partition_send[0] = partition_send[jproc].ncell[0];
				tmp_info.j2i_ncell_partition_send[1] = partition_send[jproc].ncell[1];
				tmp_info.j2i_ncell_partition_send[2] = partition_send[jproc].ncell[2];

				tmp_info.j2i_ncell_partition_recv[0] = partition_recv[iproc].ncell[0];
				tmp_info.j2i_ncell_partition_recv[1] = partition_recv[iproc].ncell[1];
				tmp_info.j2i_ncell_partition_recv[2] = partition_recv[iproc].ncell[2];

/*
//				if( rank == 0 && ncell_send[0] != ncell_recv[0] )
				if( rank == 0 )
				{
					std::cout << "ERROR!! " << std::endl;

					std::cout << "proc: " << jproc << " " << iproc << std::endl;

					std::cout << "imin: " << jmin[0] << " " << imin[0] << " --> " <<
					tmp_info.j2i_min_send[0] << " " << tmp_info.j2i_min_recv[0] << std::endl;


					std::cout << "imax: " << jmax[0] << " " << imax[0] << " --> " <<
					tmp_info.j2i_max_send[0] << " " << tmp_info.j2i_max_recv[0] << std::endl;

					std::cout << "ncells: " << ncell_send[0] << " " << ncell_recv[0] << std::endl;
				}
*/

			}


//----------------------------------------------------------------------------//
// Store communication
//----------------------------------------------------------------------------//
			if(tmp_info.nway > 0)
			{
				if( iproc == jproc ){ tmp_info.nway = 0; }
				++icomm;
				tmp_info_all.push_back( tmp_info );
			}


		} //jproc
	} //iproc

// Store total number of communications
	ncomm = icomm;

//----------------------------------------------------------------------------//
// Find communication sequence (a greedy colouring algorithm)
//----------------------------------------------------------------------------//
	proc_busy.resize(nproc, 0);
	comm_list.resize(ncomm, 1);

	class_communication::ncomm = 0;

	sum_comm = ncomm;
	while( sum_comm > 0 )
	{

/*
		if(rank == 0)
		{
			std::cout <<  "============" << std::endl;
		}
*/

// Free processors
		for( i = 0; i < nproc; ++i)
		{
			proc_busy[i] = 0;
		}

// Loop through all communications
		sum_comm = 0;
		for( i = 0; i < ncomm; ++i)
		{

			iproc = tmp_info_all[i].iproc;
			jproc = tmp_info_all[i].jproc;

// Check if communication is still listed and if processors are available
			if( comm_list[i] == 1 && proc_busy[iproc] == 0 && proc_busy[jproc] == 0 )
			{
				comm_list[i] = 0;
				proc_busy[iproc] = 1;
				proc_busy[jproc] = 1;

/*
	if(rank == 0)
	{
		std::cout << "Procs: " << iproc << " " << jproc 
		          << " - " << tmp_info_all[i].nway << std::endl;
	}
*/

// Append communication to the processor specific communication list
				if( rank == iproc || rank == jproc )
				{
					class_communication::info.push_back( tmp_info_all[i] );
					++class_communication::ncomm;
				}

			}

// Sum remaining communications
			sum_comm = sum_comm + comm_list[i];

		}

	}

//----------------------------------------------------------------------------//
// Output
//----------------------------------------------------------------------------//
/*
for (j=0; j<nproc; ++j)
{
	if(rank == j)
	{
		std::cout << "rank: " << rank << std::endl;
		std::cout << "total ncoms: " << ncomm  << std::endl;
		std::cout << "ncom rounds: " << class_communication::ncomm  << std::endl;
i = 0;
		for(i = 0; i < class_communication::ncomm; i++ )
		{
			std::cout << "imin: " << class_communication::info[i].i2j_min_send[0] 
			               << " " << class_communication::info[i].i2j_min_recv[0]
			               << " " << class_communication::info[i].i2j_min_send[1]
			               << " " << class_communication::info[i].i2j_min_recv[1]
			               << " " << class_communication::info[i].i2j_min_send[2]
			               << " " << class_communication::info[i].i2j_min_recv[2] 
			               << std::endl;
			std::cout << "imax: " << class_communication::info[i].i2j_max_send[0] 
			               << " " << class_communication::info[i].i2j_max_recv[0]
			               << " " << class_communication::info[i].i2j_max_send[1]
			               << " " << class_communication::info[i].i2j_max_recv[1]
			               << " " << class_communication::info[i].i2j_max_send[2]
			               << " " << class_communication::info[i].i2j_max_recv[2] 
			               << std::endl;
			std::cout << "--" << std::endl;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
*/
	return;
}
