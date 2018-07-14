//----------------------------------------------------------------------------//
/*
  File:         communication_setup2d.cpp

  Description:  
*/
//----------------------------------------------------------------------------//

void class_communication::setup2d( std::vector<class_partition_info> partition_send, 
                                   std::vector<class_partition_info> partition_recv )
{

//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int i,j;
	int nproc, rank;
	int iproc, jproc;
	int icomm, ncomm, sum_comm;
	int imin[2], imax[2], jmin[2], jmax[2];
	int ncell_send[2], ncell_recv[2];

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

		class_communication::partition_send.push_back( prt_info );

		prt_info.ncell[0] = partition_recv[i].ncell[0];
		prt_info.ncell[1] = partition_recv[i].ncell[1];

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

			imax[0] = partition_send[iproc].icell[0] + partition_send[iproc].ncell[0] - 1;
			imax[1] = partition_send[iproc].icell[1] + partition_send[iproc].ncell[1] - 1;

			jmin[0] = partition_recv[jproc].icell[0];
			jmin[1] = partition_recv[jproc].icell[1];

			jmax[0] = partition_recv[jproc].icell[0] + partition_recv[jproc].ncell[0] - 1;
			jmax[1] = partition_recv[jproc].icell[1] + partition_recv[jproc].ncell[1] - 1;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Check overlap
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			if( jmax[0] >= imin[0] && jmax[1] >= imin[1] &&
			    imax[0] >= jmin[0] && imax[1] >= jmin[1] )
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

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store number of cells to be communicated
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				ncell_send[0] = tmp_info.i2j_max_send[0] - tmp_info.i2j_min_send[0] + 1;
				ncell_send[1] = tmp_info.i2j_max_send[1] - tmp_info.i2j_min_send[1] + 1;

				ncell_recv[0] = tmp_info.i2j_max_recv[0] - tmp_info.i2j_min_recv[0] + 1;
				ncell_recv[1] = tmp_info.i2j_max_recv[1] - tmp_info.i2j_min_recv[1] + 1;

				if( ncell_send[0] != ncell_recv[0] || ncell_send[1] != ncell_recv[1] )
				{
					std::cout << "ERROR!! " << std::endl;
				}

				tmp_info.i2j_ncell[0] = ncell_send[0];
				tmp_info.i2j_ncell[1] = ncell_send[1];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store total number of cells of the two partition
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				tmp_info.i2j_ncell_partition_send[0] = partition_send[iproc].ncell[0];
				tmp_info.i2j_ncell_partition_send[1] = partition_send[iproc].ncell[1];

				tmp_info.i2j_ncell_partition_recv[0] = partition_recv[jproc].ncell[0];
				tmp_info.i2j_ncell_partition_recv[1] = partition_recv[jproc].ncell[1];

			}

//----------------------------------------------------------------------------//
// Overlap of subdomains from jproc to iproc
//----------------------------------------------------------------------------//
			jmin[0] = partition_send[jproc].icell[0];
			jmin[1] = partition_send[jproc].icell[1];

			jmax[0] = partition_send[jproc].icell[0] + partition_send[jproc].ncell[0] - 1;
			jmax[1] = partition_send[jproc].icell[1] + partition_send[jproc].ncell[1] - 1;

			imin[0] = partition_recv[iproc].icell[0];
			imin[1] = partition_recv[iproc].icell[1];

			imax[0] = partition_recv[iproc].icell[0] + partition_recv[iproc].ncell[0] - 1;
			imax[1] = partition_recv[iproc].icell[1] + partition_recv[iproc].ncell[1] - 1;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Check overlap
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			if( imax[0] >= jmin[0] && imax[1] >= jmin[1] &&
			    jmax[0] >= imin[0] && jmax[1] >= imin[1] &&
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

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store number of cells to be communicated
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				ncell_send[0] = tmp_info.j2i_max_send[0] - tmp_info.j2i_min_send[0] + 1;
				ncell_send[1] = tmp_info.j2i_max_send[1] - tmp_info.j2i_min_send[1] + 1;

				ncell_recv[0] = tmp_info.j2i_max_recv[0] - tmp_info.j2i_min_recv[0] + 1;
				ncell_recv[1] = tmp_info.j2i_max_recv[1] - tmp_info.j2i_min_recv[1] + 1;

				if( ncell_send[0] != ncell_recv[0] || ncell_send[1] != ncell_recv[1] )
				{
					std::cout << "ERROR!! " << std::endl;
				}

				tmp_info.j2i_ncell[0] = ncell_send[0];
				tmp_info.j2i_ncell[1] = ncell_send[1];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store total number of cells of the two partition
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				tmp_info.j2i_ncell_partition_send[0] = partition_send[jproc].ncell[0];
				tmp_info.j2i_ncell_partition_send[1] = partition_send[jproc].ncell[1];

				tmp_info.j2i_ncell_partition_recv[0] = partition_recv[iproc].ncell[0];
				tmp_info.j2i_ncell_partition_recv[1] = partition_recv[iproc].ncell[1];

/*
				if( rank == 0 && ncell_send[0] != ncell_recv[0] )
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
	if(rank == 0)
	{
		std::cout << ""  << std::endl;
		std::cout << "total ncoms: " << ncomm  << std::endl;
		std::cout << "ncom rounds: " << class_communication::ncomm  << std::endl;
		std::cout << ""  << std::endl;
	}
*/


	return;
}