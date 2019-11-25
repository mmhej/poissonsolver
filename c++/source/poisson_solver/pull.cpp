//----------------------------------------------------------------------------//
/*
  File:         pull.cpp

  Description:  Pulls the field arrays from the poisson solver object
*/
//----------------------------------------------------------------------------//

void poisson_solver::pull( double * real_rhsX, double * real_rhsY,
                           double * real_rhsZ, double * real_lhsX,
                           double * real_lhsY, double * real_lhsZ )
{

//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int i, j, k;
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

//----------------------------------------------------------------------------//
// Objects
//----------------------------------------------------------------------------//
	double *  buffer_send = NULL;
	double ** buffer_recv = new double * [xpen2real.ncomm];

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//----------------------------------------------------------------------------//
// Find out which fields to map
//----------------------------------------------------------------------------//
	nmap = 0;
	if( real_rhsX != NULL && rhsX != NULL ){ rX = true; nmap += 1; }
	if( real_rhsY != NULL && rhsY != NULL ){ rY = true; nmap += 1; }
	if( real_rhsZ != NULL && rhsZ != NULL ){ rZ = true; nmap += 1; }

	if( real_lhsX != NULL && lhsX != NULL ){ lX = true; nmap += 1; }
	if( real_lhsY != NULL && lhsY != NULL ){ lY = true; nmap += 1; }
	if( real_lhsZ != NULL && lhsZ != NULL ){ lZ = true; nmap += 1; }

//----------------------------------------------------------------------------//
// Allocate and pack communication buffers
//----------------------------------------------------------------------------//
	for( i = 0; i < xpen2real.ncomm; ++i )
	{

//----------------------------------------------------------------------------//
// Self-communication
//----------------------------------------------------------------------------//
		if( xpen2real.info[i].nway == 0 )
		{

			imin[0] = xpen2real.info[i].i2j_min_send[0];
			imin[1] = xpen2real.info[i].i2j_min_send[1];
			imin[2] = xpen2real.info[i].i2j_min_send[2];

			imax[0] = xpen2real.info[i].i2j_max_send[0];
			imax[1] = xpen2real.info[i].i2j_max_send[1];
			imax[2] = xpen2real.info[i].i2j_max_send[2];

			nx = xpen2real.info[i].i2j_ncell_partition_send[0];
			ny = xpen2real.info[i].i2j_ncell_partition_send[1];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Fill receive buffer directly
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			nsend = nmap * xpen2real.info[i].i2j_ncell[0]
			             * xpen2real.info[i].i2j_ncell[1]
			             * xpen2real.info[i].i2j_ncell[2];
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
							buffer_recv[i][k] = std::real(rhsX[pqs]);
							++k;
						}
						if( rY )
						{
							buffer_recv[i][k] = std::real(rhsY[pqs]);
							++k;
						}
						if( rZ )
						{
							buffer_recv[i][k] = std::real(rhsZ[pqs]);
							++k;
						}

						if( lX )
						{
							buffer_recv[i][k] = std::real(lhsX[pqs]);
							++k;
						}
						if( lY )
						{
							buffer_recv[i][k] = std::real(lhsY[pqs]);
							++k;
						}
						if( lZ )
						{
							buffer_recv[i][k] = std::real(lhsZ[pqs]);
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

			iproc = xpen2real.info[i].iproc;
			jproc = xpen2real.info[i].jproc;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Mesh info
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			pack = true;
			if( rank == iproc && xpen2real.info[i].nway != 2 )
			{

				imin[0] = xpen2real.info[i].i2j_min_send[0];
				imin[1] = xpen2real.info[i].i2j_min_send[1];
				imin[2] = xpen2real.info[i].i2j_min_send[2];

				imax[0] = xpen2real.info[i].i2j_max_send[0];
				imax[1] = xpen2real.info[i].i2j_max_send[1];
				imax[2] = xpen2real.info[i].i2j_max_send[2];

				nsend = nmap * xpen2real.info[i].i2j_ncell[0]
				             * xpen2real.info[i].i2j_ncell[1]
				             * xpen2real.info[i].i2j_ncell[2];

				nx = xpen2real.info[i].i2j_ncell_partition_send[0];
				ny = xpen2real.info[i].i2j_ncell_partition_send[1];
			}
			else if( rank == jproc && xpen2real.info[i].nway > 1 )
			{

				imin[0] = xpen2real.info[i].j2i_min_send[0];
				imin[1] = xpen2real.info[i].j2i_min_send[1];
				imin[2] = xpen2real.info[i].j2i_min_send[2];

				imax[0] = xpen2real.info[i].j2i_max_send[0];
				imax[1] = xpen2real.info[i].j2i_max_send[1];
				imax[2] = xpen2real.info[i].j2i_max_send[2];

				nsend = nmap * xpen2real.info[i].j2i_ncell[0]
				             * xpen2real.info[i].j2i_ncell[1]
				             * xpen2real.info[i].j2i_ncell[2];

				nx = xpen2real.info[i].j2i_ncell_partition_send[0];
				ny = xpen2real.info[i].j2i_ncell_partition_send[1];
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
								buffer_send[k] = std::real(rhsX[pqs]);
								++k;
							}
							if( rY )
							{
								buffer_send[k] = std::real(rhsY[pqs]);
								++k;
							}
							if( rZ )
							{
								buffer_send[k] = std::real(rhsZ[pqs]);
								++k;
							}
							if( lX )
							{
								buffer_send[k] = std::real(lhsX[pqs]);
								++k;
							}
							if( lY )
							{
								buffer_send[k] = std::real(lhsY[pqs]);
								++k;
							}
							if( lZ )
							{
								buffer_send[k] = std::real(lhsZ[pqs]);
								++k;
							}

						}
					}
				}

			}

//----------------------------------------------------------------------------//
// Send/Recieve
//----------------------------------------------------------------------------//
			if( xpen2real.info[i].nway == 1 ) // only i to j
			{
				if( rank == iproc ) // Sender
				{
					MPI_Send( buffer_send, nsend, MPI_DOUBLE, jproc, jproc,
					          MPI_COMM_WORLD );
				}
				else if( rank == jproc ) // Receiver
				{
					nrecv = nmap * xpen2real.info[i].i2j_ncell[0] 
					             * xpen2real.info[i].i2j_ncell[1] 
					             * xpen2real.info[i].i2j_ncell[2];
					buffer_recv[i] = new double[ nrecv ];
					MPI_Recv( buffer_recv[i], nrecv, MPI_DOUBLE, iproc, jproc, 
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
			}
			else if( xpen2real.info[i].nway == 2 ) // only j to i
			{
				if( rank == jproc ) // Sender
				{
					MPI_Send( buffer_send, nsend, MPI_DOUBLE, iproc, iproc,
					          MPI_COMM_WORLD );
				}
				else if( rank == iproc ) // Receiver
				{
					nrecv = nmap * xpen2real.info[i].j2i_ncell[0] 
					             * xpen2real.info[i].j2i_ncell[1] 
					             * xpen2real.info[i].j2i_ncell[2];
					buffer_recv[i] = new double[ nrecv ];
					MPI_Recv( buffer_recv[i], nrecv, MPI_DOUBLE, jproc, iproc, 
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
			}
			else if( xpen2real.info[i].nway == 3 ) // i to j and j to i
			{
				if( rank == iproc )
				{
					nrecv = nmap * xpen2real.info[i].j2i_ncell[0]
					             * xpen2real.info[i].j2i_ncell[1]
					             * xpen2real.info[i].j2i_ncell[2];
					buffer_recv[i] = new double[ nrecv ];
					MPI_Sendrecv( buffer_send, nsend, MPI_DOUBLE, jproc, iproc, 
					              buffer_recv[i], nrecv, MPI_DOUBLE, jproc, jproc, 
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				else if( rank == jproc )
				{
					nrecv = nmap * xpen2real.info[i].i2j_ncell[0]
					             * xpen2real.info[i].i2j_ncell[1]
					             * xpen2real.info[i].i2j_ncell[2];
					buffer_recv[i] = new double[ nrecv ];
					MPI_Sendrecv( buffer_send, nsend, MPI_DOUBLE, iproc, jproc, 
					              buffer_recv[i], nrecv, MPI_DOUBLE, iproc, iproc, 
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			} // nway

		}

//----------------------------------------------------------------------------//
// De-allocate send buffers
//----------------------------------------------------------------------------//
		if( buffer_send != NULL)
		{
			delete [] buffer_send;
			buffer_send = NULL;
		}

//		MPI_Barrier(MPI_COMM_WORLD);
	} // ncomm

//----------------------------------------------------------------------------//
// Clear arrays
//----------------------------------------------------------------------------//
	nrecv = xpen2real.partition_recv[rank].ncell[0] 
	      * xpen2real.partition_recv[rank].ncell[1] 
	      * xpen2real.partition_recv[rank].ncell[2];

/*
	if( rX ){ real_rhsX = new double[nrecv]; }//else{ real_rhsX = NULL; }
	if( rY ){ real_rhsY = new double[nrecv]; }//else{ real_rhsY = NULL; }
	if( rZ ){ real_rhsZ = new double[nrecv]; }//else{ real_rhsZ = NULL; }

	if( lX ){ real_lhsX = new double[nrecv]; }//else{ real_lhsX = NULL; }
	if( lY ){ real_lhsY = new double[nrecv]; }//else{ real_lhsY = NULL; }
	if( lZ ){ real_lhsZ = new double[nrecv]; }//else{ real_lhsZ = NULL; }
*/

	for(i = 0; i < nrecv; i++)
	{
		if( rX ){ real_rhsX[i] = 0.0; }//else{ real_rhsX = NULL; }
		if( rY ){ real_rhsY[i] = 0.0; }//else{ real_rhsY = NULL; }
		if( rZ ){ real_rhsZ[i] = 0.0; }//else{ real_rhsZ = NULL; }

		if( lX ){ real_lhsX[i] = 0.0; }//else{ real_lhsX = NULL; }
		if( lY ){ real_lhsY[i] = 0.0; }//else{ real_lhsY = NULL; }
		if( lZ ){ real_lhsZ[i] = 0.0; }//else{ real_lhsZ = NULL; }
	}

//----------------------------------------------------------------------------//
// Unpack communication buffer
//----------------------------------------------------------------------------//
	for( i = 0; i < xpen2real.ncomm; ++i )
	{

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Mesh info
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		pack = true;
		iproc = xpen2real.info[i].iproc;
		jproc = xpen2real.info[i].jproc;

		if( rank == jproc && xpen2real.info[i].nway != 2 )
		{
			imin[0] = xpen2real.info[i].i2j_min_recv[0];
			imin[1] = xpen2real.info[i].i2j_min_recv[1];
			imin[2] = xpen2real.info[i].i2j_min_recv[2];

			imax[0] = xpen2real.info[i].i2j_max_recv[0];
			imax[1] = xpen2real.info[i].i2j_max_recv[1];
			imax[2] = xpen2real.info[i].i2j_max_recv[2];

			nx = xpen2real.info[i].i2j_ncell_partition_recv[0];
			ny = xpen2real.info[i].i2j_ncell_partition_recv[1];
		}
		else if( rank == iproc && xpen2real.info[i].nway > 1 )
		{
			imin[0] = xpen2real.info[i].j2i_min_recv[0];
			imin[1] = xpen2real.info[i].j2i_min_recv[1];
			imin[2] = xpen2real.info[i].j2i_min_recv[2];

			imax[0] = xpen2real.info[i].j2i_max_recv[0];
			imax[1] = xpen2real.info[i].j2i_max_recv[1];
			imax[2] = xpen2real.info[i].j2i_max_recv[2];

			nx = xpen2real.info[i].j2i_ncell_partition_recv[0];
			ny = xpen2real.info[i].j2i_ncell_partition_recv[1];
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
							real_rhsX[pqs] = buffer_recv[i][k];
							++k;
						}
						if( rY )
						{
							real_rhsY[pqs] = buffer_recv[i][k];
							++k;
						}
						if( rZ )
						{
							real_rhsZ[pqs] = buffer_recv[i][k];
							++k;
						}
						if( lX )
						{
							real_lhsX[pqs] = buffer_recv[i][k];
							++k;
						}
						if( lY )
						{
							real_lhsY[pqs] = buffer_recv[i][k];
							++k;
						}
						if( lZ )
						{
							real_lhsZ[pqs] = buffer_recv[i][k];
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
		for( i = 0; i < xpen2real.ncomm; ++i )
		{
			delete [] buffer_recv[i];
			buffer_recv[i] = NULL;
		}
		delete [] buffer_recv;
		buffer_recv = NULL;
	}

	return;
}
