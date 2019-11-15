//----------------------------------------------------------------------------//
/*
  File:         map.cpp

  Description:  Performs the MPI mapping of the input communication object
*/
//----------------------------------------------------------------------------//

void poisson_solver::map( class_communication comm )
{

//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int i,j,k;
	int p, q, s, sn, sqn, pqs;
	int nx, ny;
	int nproc, rank;
	int iproc, jproc;
	int nsend, nrecv;
	int imin[3], imax[3];
	int ncell[3];
	int nmap;
	bool pack;

	bool rX = false;
	bool rY = false;
	bool rZ = false;

	bool lX = false;
	bool lY = false;
	bool lZ = false;

	bool mG = false;

//----------------------------------------------------------------------------//
// Objects
//----------------------------------------------------------------------------//
	double *  buffer_send = NULL;
	double ** buffer_recv = new double * [comm.ncomm];

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//----------------------------------------------------------------------------//
// Find out which fields to map
//----------------------------------------------------------------------------//
	nmap = 0;
	if( rhsX != NULL ){ rX = true; nmap += 1; }
	if( rhsY != NULL ){ rY = true; nmap += 1; }
	if( rhsZ != NULL ){ rZ = true; nmap += 1; }

	if( lhsX != NULL ){ lX = true; nmap += 1; }
	if( lhsY != NULL ){ lY = true; nmap += 1; }
	if( lhsZ != NULL ){ lZ = true; nmap += 1; }

	if( mapG != NULL ){ mG = true; nmap += 1; }

//----------------------------------------------------------------------------//
// Allocate and pack communication buffers
//----------------------------------------------------------------------------//
	for( i = 0; i < comm.ncomm; ++i )
	{

//----------------------------------------------------------------------------//
// Self-communication
//----------------------------------------------------------------------------//
		if( comm.info[i].nway == 0 )
		{

			imin[0] = comm.info[i].i2j_min_send[0];
			imin[1] = comm.info[i].i2j_min_send[1];
			imin[2] = comm.info[i].i2j_min_send[2];

			imax[0] = comm.info[i].i2j_max_send[0];
			imax[1] = comm.info[i].i2j_max_send[1];
			imax[2] = comm.info[i].i2j_max_send[2];

			nx = comm.info[i].i2j_ncell_partition_send[0];
			ny = comm.info[i].i2j_ncell_partition_send[1];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Fill receive buffer directly
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			nsend = 2 * nmap * comm.info[i].i2j_ncell[0]
			                 * comm.info[i].i2j_ncell[1]
			                 * comm.info[i].i2j_ncell[2];

//std::cout << "nsend " << nsend << " " << imin[0] << " " << imax[0] << " " << imin[1] << " " << imax[1] << " " << imin[2] << " " << imax[2] << std::endl;
			buffer_recv[i] = new double[ nsend ]();

			k = 0;
			for( s = imin[2]; s <= imax[2]; ++s )
			{
				sn = s * ny;
				for( q = imin[1]; q <= imax[1]; ++q )
				{
					sqn = ( sn + q ) * nx;
					for( p = imin[0]; p <= imax[0]; ++p )
					{
						pqs = sqn + p;

						if( rX )
						{
							buffer_recv[i][k] = std::real(rhsX[pqs]);
							++k;
							buffer_recv[i][k] = std::imag(rhsX[pqs]);
							++k;
						}
						if( rY )
						{
							buffer_recv[i][k] = std::real(rhsY[pqs]);
							++k;
							buffer_recv[i][k] = std::imag(rhsY[pqs]);
							++k;
						}
						if( rZ )
						{
							buffer_recv[i][k] = std::real(rhsZ[pqs]);
							++k;
							buffer_recv[i][k] = std::imag(rhsZ[pqs]);
							++k;
						}

						if( lX )
						{
							buffer_recv[i][k] = std::real(lhsX[pqs]);
							++k;
							buffer_recv[i][k] = std::imag(lhsX[pqs]);
							++k;
						}
						if( lY )
						{
							buffer_recv[i][k] = std::real(lhsY[pqs]);
							++k;
							buffer_recv[i][k] = std::imag(lhsY[pqs]);
							++k;
						}
						if( lZ )
						{
							buffer_recv[i][k] = std::real(lhsZ[pqs]);
							++k;
							buffer_recv[i][k] = std::imag(lhsZ[pqs]);
							++k;
						}

						if( mG )
						{
							buffer_recv[i][k] = std::real(mapG[pqs]);
							++k;
							buffer_recv[i][k] = std::imag(mapG[pqs]);
							++k;
						}

					}
				}
			}

		}

//----------------------------------------------------------------------------//
// Inter-communication
//----------------------------------------------------------------------------//
		else
		{

			iproc = comm.info[i].iproc;
			jproc = comm.info[i].jproc;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Mesh info
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			pack = true;
			if( rank == iproc && comm.info[i].nway != 2 )
			{

				imin[0] = comm.info[i].i2j_min_send[0];
				imin[1] = comm.info[i].i2j_min_send[1];
				imin[2] = comm.info[i].i2j_min_send[2];

				imax[0] = comm.info[i].i2j_max_send[0];
				imax[1] = comm.info[i].i2j_max_send[1];
				imax[2] = comm.info[i].i2j_max_send[2];

				nsend = 2 * nmap * comm.info[i].i2j_ncell[0]
				                 * comm.info[i].i2j_ncell[1]
				                 * comm.info[i].i2j_ncell[2];

				nx = comm.info[i].i2j_ncell_partition_send[0];
				ny = comm.info[i].i2j_ncell_partition_send[1];
			}
			else if( rank == jproc && comm.info[i].nway > 1 )
			{

				imin[0] = comm.info[i].j2i_min_send[0];
				imin[1] = comm.info[i].j2i_min_send[1];
				imin[2] = comm.info[i].j2i_min_send[2];

				imax[0] = comm.info[i].j2i_max_send[0];
				imax[1] = comm.info[i].j2i_max_send[1];
				imax[2] = comm.info[i].j2i_max_send[2];

				nsend = 2 * nmap * comm.info[i].j2i_ncell[0]
				                 * comm.info[i].j2i_ncell[1]
				                 * comm.info[i].j2i_ncell[2];

				nx = comm.info[i].j2i_ncell_partition_send[0];
				ny = comm.info[i].j2i_ncell_partition_send[1];
			}
			else
			{
				pack = false;
			}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Pack send buffer
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			if(pack)
			{
				buffer_send = new double[ nsend ]();
				k = 0;

				for( s = imin[2]; s <= imax[2]; ++s )
				{
					sn = s * ny;
					for( q = imin[1]; q <= imax[1]; ++q )
					{
						sqn = ( sn + q ) * nx;
						for( p = imin[0]; p <= imax[0]; ++p )
						{
							pqs = sqn + p;


							if( rX )
							{
								buffer_send[k] = std::real(rhsX[pqs]);
								++k;
								buffer_send[k] = std::imag(rhsX[pqs]);
								++k;
							}
							if( rY )
							{
								buffer_send[k] = std::real(rhsY[pqs]);
								++k;
								buffer_send[k] = std::imag(rhsY[pqs]);
								++k;
							}
							if( rZ )
							{
								buffer_send[k] = std::real(rhsZ[pqs]);
								++k;
								buffer_send[k] = std::imag(rhsZ[pqs]);
								++k;
							}
							if( lX )
							{
								buffer_send[k] = std::real(lhsX[pqs]);
								++k;
								buffer_send[k] = std::imag(lhsX[pqs]);
								++k;
							}
							if( lY )
							{
								buffer_send[k] = std::real(lhsY[pqs]);
								++k;
								buffer_send[k] = std::imag(lhsY[pqs]);
								++k;
							}
							if( lZ )
							{
								buffer_send[k] = std::real(lhsZ[pqs]);
								++k;
								buffer_send[k] = std::imag(lhsZ[pqs]);
								++k;
							}

							if( mG )
							{
								buffer_send[k] = std::real(mapG[pqs]);
								++k;
								buffer_send[k] = std::imag(mapG[pqs]);
								++k;
							}

						}
					}
				}

			}

//----------------------------------------------------------------------------//
// Send/Recieve
//----------------------------------------------------------------------------//
			if( comm.info[i].nway == 1 ) // only i to j
			{
				if( rank == iproc ) // Sender
				{
					MPI_Send( buffer_send, nsend, MPI_DOUBLE, jproc, jproc,
					          MPI_COMM_WORLD );
				}
				else if( rank == jproc ) // Receiver
				{
					nrecv = 2 * nmap * comm.info[i].i2j_ncell[0] 
					                 * comm.info[i].i2j_ncell[1] 
					                 * comm.info[i].i2j_ncell[2];
					buffer_recv[i] = new double[ nrecv ]();
					MPI_Recv( buffer_recv[i], nrecv, MPI_DOUBLE, iproc, jproc, 
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
			}
			else if( comm.info[i].nway == 2 ) // only j to i
			{
				if( rank == jproc ) // Sender
				{
					MPI_Send( buffer_send, nsend, MPI_DOUBLE, iproc, iproc,
					          MPI_COMM_WORLD );
				}
				else if( rank == iproc ) // Receiver
				{
					nrecv = 2 * nmap * comm.info[i].j2i_ncell[0] 
					                 * comm.info[i].j2i_ncell[1] 
					                 * comm.info[i].j2i_ncell[2];
					buffer_recv[i] = new double[ nrecv ]();
					MPI_Recv( buffer_recv[i], nrecv, MPI_DOUBLE, jproc, iproc, 
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
			}
			else if( comm.info[i].nway == 3 ) // i to j and j to i
			{
				if( rank == iproc )
				{
					nrecv = 2 * nmap * comm.info[i].j2i_ncell[0]
					                 * comm.info[i].j2i_ncell[1]
					                 * comm.info[i].j2i_ncell[2];
					buffer_recv[i] = new double[ nrecv ]();
					MPI_Sendrecv( buffer_send, nsend, MPI_DOUBLE, jproc, iproc, 
					              buffer_recv[i], nrecv, MPI_DOUBLE, jproc, jproc, 
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				else if( rank == jproc )
				{
					nrecv = 2 * nmap * comm.info[i].i2j_ncell[0]
					                 * comm.info[i].i2j_ncell[1]
					                 * comm.info[i].i2j_ncell[2];
					buffer_recv[i] = new double[ nrecv ]();
					MPI_Sendrecv( buffer_send, nsend, MPI_DOUBLE, iproc, jproc, 
					              buffer_recv[i], nrecv, MPI_DOUBLE, iproc, iproc, 
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			} // nway

		}

//		MPI_Barrier(MPI_COMM_WORLD);
	} // ncomm

//----------------------------------------------------------------------------//
// De-allocate send buffers
//----------------------------------------------------------------------------//
	if( buffer_send != NULL)
	{
		delete [] buffer_send;
		buffer_send = NULL;
	}

//----------------------------------------------------------------------------//
// Resize arrays
//----------------------------------------------------------------------------//
	nrecv = comm.partition_recv[rank].ncell[0] 
	      * comm.partition_recv[rank].ncell[1] 
	      * comm.partition_recv[rank].ncell[2];

	if( rX ){ rhsX = new std::complex<double>[nrecv](); }
	if( rY ){ rhsY = new std::complex<double>[nrecv](); }
	if( rZ ){ rhsZ = new std::complex<double>[nrecv](); }

	if( lX ){ lhsX = new std::complex<double>[nrecv](); }
	if( lY ){ lhsY = new std::complex<double>[nrecv](); }
	if( lZ ){ lhsZ = new std::complex<double>[nrecv](); }

	if( mG ){ mapG = new std::complex<double>[nrecv](); }

//----------------------------------------------------------------------------//
// Unpack communication buffer
//----------------------------------------------------------------------------//
	for( i = 0; i < comm.ncomm; ++i )
	{

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Mesh info
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		pack = true;
		iproc = comm.info[i].iproc;
		jproc = comm.info[i].jproc;

		if( rank == jproc && comm.info[i].nway != 2 )
		{
			imin[0] = comm.info[i].i2j_min_recv[0];
			imin[1] = comm.info[i].i2j_min_recv[1];
			imin[2] = comm.info[i].i2j_min_recv[2];

			imax[0] = comm.info[i].i2j_max_recv[0];
			imax[1] = comm.info[i].i2j_max_recv[1];
			imax[2] = comm.info[i].i2j_max_recv[2];

			nx = comm.info[i].i2j_ncell_partition_recv[0];
			ny = comm.info[i].i2j_ncell_partition_recv[1];
		}
		else if( rank == iproc && comm.info[i].nway > 1 )
		{
			imin[0] = comm.info[i].j2i_min_recv[0];
			imin[1] = comm.info[i].j2i_min_recv[1];
			imin[2] = comm.info[i].j2i_min_recv[2];

			imax[0] = comm.info[i].j2i_max_recv[0];
			imax[1] = comm.info[i].j2i_max_recv[1];
			imax[2] = comm.info[i].j2i_max_recv[2];

			nx = comm.info[i].j2i_ncell_partition_recv[0];
			ny = comm.info[i].j2i_ncell_partition_recv[1];
		}
		else
		{
			pack = false;
		}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Unpack buffers
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		if(pack)
		{
			k = 0;

			for( s = imin[2]; s <= imax[2]; ++s )
			{
				sn = s * ny;
				for( q = imin[1]; q <= imax[1]; ++q )
				{
					sqn = ( sn + q ) * nx;
					for( p = imin[0]; p <= imax[0]; ++p )
					{
						pqs = sqn + p;

						if( rX )
						{
							rhsX[pqs] = { buffer_recv[i][k], buffer_recv[i][k+1] };
							k += 2;
						}
						if( rY )
						{
							rhsY[pqs] = { buffer_recv[i][k], buffer_recv[i][k+1] };
							k += 2;
						}
						if( rZ )
						{
							rhsZ[pqs] = { buffer_recv[i][k], buffer_recv[i][k+1] };
							k += 2;
						}
						if( lX )
						{
							lhsX[pqs] = { buffer_recv[i][k], buffer_recv[i][k+1] };
							k += 2;
						}
						if( lY )
						{
							lhsY[pqs] = { buffer_recv[i][k], buffer_recv[i][k+1] };
							k += 2;
						}
						if( lZ )
						{
							lhsZ[pqs] = { buffer_recv[i][k], buffer_recv[i][k+1] };
							k += 2;
						}

						if( mG )
						{
							mapG[pqs] = { buffer_recv[i][k], buffer_recv[i][k+1] };
							k += 2;
						}

					}
				}
			}

		}

	} // ncomm


//----------------------------------------------------------------------------//
// De-allocate recieve buffers
//----------------------------------------------------------------------------//
	if( buffer_recv != NULL)
	{
		for( i = 0; i < comm.ncomm; ++i )
		{
			delete [] buffer_recv[i];
			buffer_recv[i] = NULL;
		}
		delete [] buffer_recv;
		buffer_recv = NULL;
	}

	return;
}
