//----------------------------------------------------------------------------//
/*
  File:         push.cpp

  Description:  
*/
//----------------------------------------------------------------------------//

void class_greenfish::push( double * real_rhsX, double * real_rhsY,
                            double * real_rhsZ, double * real_lhsX,
                            double * real_lhsY, double * real_lhsZ )
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
	int nmap;
	bool pack;

	bool rX = false;
	bool rY = false;
	bool rZ = false;
	bool lX = false;
	bool lY = false;
	bool lZ = false;

//----------------------------------------------------------------------------//
// Objects
//----------------------------------------------------------------------------//
	double *  buffer_send = NULL;
	double ** buffer_recv = new double * [real2xpen.ncomm];

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//----------------------------------------------------------------------------//
// Find out which fields to map
//----------------------------------------------------------------------------//
	nmap = 0;
	if( real_rhsX != NULL ){ rX = true; nmap += 1; }
	if( real_rhsY != NULL ){ rY = true; nmap += 1; }
	if( real_rhsZ != NULL ){ rZ = true; nmap += 1; }

	if( real_lhsX != NULL ){ lX = true; nmap += 1; }
	if( real_lhsY != NULL ){ lY = true; nmap += 1; }
	if( real_lhsZ != NULL ){ lZ = true; nmap += 1; }

//----------------------------------------------------------------------------//
// Allocate and pack communication buffers
//----------------------------------------------------------------------------//
	for( i = 0; i < real2xpen.ncomm; ++i )
	{

//----------------------------------------------------------------------------//
// Self-communication
//----------------------------------------------------------------------------//
		if( real2xpen.info[i].nway == 0 )
		{

			imin[0] = real2xpen.info[i].i2j_min_send[0];
			imin[1] = real2xpen.info[i].i2j_min_send[1];
			imin[2] = real2xpen.info[i].i2j_min_send[2];

			imax[0] = real2xpen.info[i].i2j_max_send[0];
			imax[1] = real2xpen.info[i].i2j_max_send[1];
			imax[2] = real2xpen.info[i].i2j_max_send[2];

			nx = real2xpen.info[i].i2j_ncell_partition_send[0];
			ny = real2xpen.info[i].i2j_ncell_partition_send[1];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Fill receive buffer directly
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			nsend = nmap * real2xpen.info[i].i2j_ncell[0]
			             * real2xpen.info[i].i2j_ncell[1]
			             * real2xpen.info[i].i2j_ncell[2];
			buffer_recv[i] = new double[ nsend ];

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
							buffer_recv[i][k] = real_rhsX[pqs];
							++k;
						}
						if( rY )
						{
							buffer_recv[i][k] = real_rhsY[pqs];
							++k;
						}
						if( rZ )
						{
							buffer_recv[i][k] = real_rhsZ[pqs];
							++k;
						}

						if( lX )
						{
							buffer_recv[i][k] = real_lhsX[pqs];
							++k;
						}
						if( lY )
						{
							buffer_recv[i][k] = real_lhsY[pqs];
							++k;
						}
						if( lZ )
						{
							buffer_recv[i][k] = real_lhsZ[pqs];
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

			iproc = real2xpen.info[i].iproc;
			jproc = real2xpen.info[i].jproc;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Mesh info
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			pack = true;
			if( rank == iproc && real2xpen.info[i].nway != 2 )
			{
				imin[0] = real2xpen.info[i].i2j_min_send[0];
				imin[1] = real2xpen.info[i].i2j_min_send[1];
				imin[2] = real2xpen.info[i].i2j_min_send[2];

				imax[0] = real2xpen.info[i].i2j_max_send[0];
				imax[1] = real2xpen.info[i].i2j_max_send[1];
				imax[2] = real2xpen.info[i].i2j_max_send[2];

				nsend = nmap * real2xpen.info[i].i2j_ncell[0]
				             * real2xpen.info[i].i2j_ncell[1]
				             * real2xpen.info[i].i2j_ncell[2];

				nx = real2xpen.info[i].i2j_ncell_partition_send[0];
				ny = real2xpen.info[i].i2j_ncell_partition_send[1];
			}
			else if( rank == jproc && real2xpen.info[i].nway > 1 )
			{
				imin[0] = real2xpen.info[i].j2i_min_send[0];
				imin[1] = real2xpen.info[i].j2i_min_send[1];
				imin[2] = real2xpen.info[i].j2i_min_send[2];

				imax[0] = real2xpen.info[i].j2i_max_send[0];
				imax[1] = real2xpen.info[i].j2i_max_send[1];
				imax[2] = real2xpen.info[i].j2i_max_send[2];

				nsend = nmap * real2xpen.info[i].j2i_ncell[0]
				             * real2xpen.info[i].j2i_ncell[1]
				             * real2xpen.info[i].j2i_ncell[2];

				nx = real2xpen.info[i].j2i_ncell_partition_send[0];
				ny = real2xpen.info[i].j2i_ncell_partition_send[1];
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
				buffer_send = new double[ nsend ];
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
								buffer_send[k] = real_rhsX[pqs];
								++k;
							}
							if( rY )
							{
								buffer_send[k] = real_rhsY[pqs];
								++k;
							}
							if( rZ )
							{
								buffer_send[k] = real_rhsZ[pqs];
								++k;
							}
							if( lX )
							{
								buffer_send[k] = real_lhsX[pqs];
								++k;
							}
							if( lY )
							{
								buffer_send[k] = real_lhsY[pqs];
								++k;
							}
							if( lZ )
							{
								buffer_send[k] = real_lhsZ[pqs];
								++k;
							}

						}
					}
				}

			}

//----------------------------------------------------------------------------//
// Send/Recieve
//----------------------------------------------------------------------------//
			if( real2xpen.info[i].nway == 1 ) // only i to j
			{
				if( rank == iproc ) // Sender
				{
					MPI_Send( buffer_send, nsend, MPI_DOUBLE, jproc, jproc,
					          MPI_COMM_WORLD );
				}
				else if( rank == jproc ) // Receiver
				{
					nrecv = nmap * real2xpen.info[i].i2j_ncell[0]
					             * real2xpen.info[i].i2j_ncell[1]
					             * real2xpen.info[i].i2j_ncell[2];
					buffer_recv[i] = new double[ nrecv ];
					MPI_Recv( buffer_recv[i], nrecv, MPI_DOUBLE, iproc, jproc, 
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
			}
			else if( real2xpen.info[i].nway == 2 ) // only j to i
			{
				if( rank == jproc ) // Sender
				{
					MPI_Send( buffer_send, nsend, MPI_DOUBLE, iproc, iproc,
					          MPI_COMM_WORLD );
				}
				else if( rank == iproc ) // Receiver
				{
					nrecv = nmap * real2xpen.info[i].j2i_ncell[0]
					             * real2xpen.info[i].j2i_ncell[1]
					             * real2xpen.info[i].j2i_ncell[2];
					buffer_recv[i] = new double[ nrecv ];
					MPI_Recv( buffer_recv[i], nrecv, MPI_DOUBLE, jproc, iproc, 
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
			}
			else if( real2xpen.info[i].nway == 3 ) // i to j and j to i
			{
				if( rank == iproc )
				{
					nrecv = nmap * real2xpen.info[i].j2i_ncell[0]
					             * real2xpen.info[i].j2i_ncell[1]
					             * real2xpen.info[i].j2i_ncell[2];
					buffer_recv[i] = new double[ nrecv ];
					MPI_Sendrecv( buffer_send, nsend, MPI_DOUBLE, jproc, iproc, 
					              buffer_recv[i], nrecv, MPI_DOUBLE, jproc, jproc, 
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				else if( rank == jproc )
				{
					nrecv = nmap * real2xpen.info[i].i2j_ncell[0]
					             * real2xpen.info[i].i2j_ncell[1]
					             * real2xpen.info[i].i2j_ncell[2];
					buffer_recv[i] = new double[ nrecv ];
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
	nrecv = real2xpen.partition_recv[rank].ncell[0]
	      * real2xpen.partition_recv[rank].ncell[1]
	      * real2xpen.partition_recv[rank].ncell[2];

	if( rX ){ rhsX = new std::complex<double>[nrecv](); }
	if( rY ){ rhsY = new std::complex<double>[nrecv](); }
	if( rZ ){ rhsZ = new std::complex<double>[nrecv](); }

	if( lX ){ lhsX = new std::complex<double>[nrecv](); }
	if( lY ){ lhsY = new std::complex<double>[nrecv](); }
	if( lZ ){ lhsZ = new std::complex<double>[nrecv](); }

//----------------------------------------------------------------------------//
// Unpack communication buffer
//----------------------------------------------------------------------------//
	for( i = 0; i < real2xpen.ncomm; ++i )
	{

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Mesh info
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		pack = true;
		iproc = real2xpen.info[i].iproc;
		jproc = real2xpen.info[i].jproc;

		if( rank == jproc && real2xpen.info[i].nway != 2 )
		{
			imin[0] = real2xpen.info[i].i2j_min_recv[0];
			imin[1] = real2xpen.info[i].i2j_min_recv[1];
			imin[2] = real2xpen.info[i].i2j_min_recv[2];

			imax[0] = real2xpen.info[i].i2j_max_recv[0];
			imax[1] = real2xpen.info[i].i2j_max_recv[1];
			imax[2] = real2xpen.info[i].i2j_max_recv[2];

			nx = real2xpen.info[i].i2j_ncell_partition_recv[0];
			ny = real2xpen.info[i].i2j_ncell_partition_recv[1];

		}
		else if( rank == iproc && real2xpen.info[i].nway > 1 )
		{
			imin[0] = real2xpen.info[i].j2i_min_recv[0];
			imin[1] = real2xpen.info[i].j2i_min_recv[1];
			imin[2] = real2xpen.info[i].j2i_min_recv[2];

			imax[0] = real2xpen.info[i].j2i_max_recv[0];
			imax[1] = real2xpen.info[i].j2i_max_recv[1];
			imax[2] = real2xpen.info[i].j2i_max_recv[2];

			nx = real2xpen.info[i].j2i_ncell_partition_recv[0];
			ny = real2xpen.info[i].j2i_ncell_partition_recv[1];
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
							rhsX[pqs] = { buffer_recv[i][k], 0.0 };
							++k;
						}
						if( rY )
						{
							rhsY[pqs] = { buffer_recv[i][k], 0.0 };
							++k;
						}
						if( rZ )
						{
							rhsZ[pqs] = { buffer_recv[i][k], 0.0 };
							++k;
						}
						if( lX )
						{
							lhsX[pqs] = { buffer_recv[i][k], 0.0 };
							++k;
						}
						if( lY )
						{
							lhsY[pqs] = { buffer_recv[i][k], 0.0 };
							++k;
						}
						if( lZ )
						{
							lhsZ[pqs] = { buffer_recv[i][k], 0.0 };
							++k;
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
		for( i = 0; i < real2xpen.ncomm; ++i )
		{
			delete [] buffer_recv[i];
			buffer_recv[i] = NULL;
		}
		delete [] buffer_recv;
		buffer_recv = NULL;
	}

//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}
