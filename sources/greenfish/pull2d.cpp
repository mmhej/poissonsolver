//----------------------------------------------------------------------------//
/*
  File:         pull2d.cpp

  Description:  
*/
//----------------------------------------------------------------------------//

void class_greenfish::pull2d( double * real_rhsX, double * real_rhsY,
                              double * real_rhsZ, double * real_lhsX,
                              double * real_lhsY, double * real_lhsZ )
{

//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int i,j,k,p,q,pq;
	int nx;
	int nproc, rank;
	int iproc, jproc;
	int nsend, nrecv;
	int imin[2], imax[2];
	int ncell[2];
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
	double * buffer_send = NULL;
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

//		if(rank == 0){ std::cout << "COMM: " << i << std::endl; }
//		MPI_Barrier(MPI_COMM_WORLD);
//----------------------------------------------------------------------------//
// Self-communication
//----------------------------------------------------------------------------//
		if( xpen2real.info[i].nway == 0 )
		{

			imin[0] = xpen2real.info[i].i2j_min_send[0];
			imin[1] = xpen2real.info[i].i2j_min_send[1];

			imax[0] = xpen2real.info[i].i2j_max_send[0];
			imax[1] = xpen2real.info[i].i2j_max_send[1];

			nx = xpen2real.info[i].i2j_ncell_partition_send[0];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Fill receive buffer directly
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			nsend = nmap * xpen2real.info[i].i2j_ncell[0]
			             * xpen2real.info[i].i2j_ncell[1];
			buffer_recv[i] = new double[ nsend ];

			k = 0;
			for( q = imin[1]; q <= imax[1]; ++q )
			{
				for( p = imin[0]; p <= imax[0]; ++p )
				{

					if( rX )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(rhsX[pq]);
						++k;
					}
					if( rY )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(rhsY[pq]);
						++k;
					}
					if( rZ )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(rhsZ[pq]);
						++k;
					}

					if( lX )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(lhsX[pq]);
						++k;
					}
					if( lY )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(lhsY[pq]);
						++k;
					}
					if( lZ )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(lhsZ[pq]);
						++k;
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

				imax[0] = xpen2real.info[i].i2j_max_send[0];
				imax[1] = xpen2real.info[i].i2j_max_send[1];

				nsend = nmap * xpen2real.info[i].i2j_ncell[0]
				             * xpen2real.info[i].i2j_ncell[1];
				nx = xpen2real.info[i].i2j_ncell_partition_send[0];
			}
			else if( rank == jproc && xpen2real.info[i].nway > 1 )
			{

				imin[0] = xpen2real.info[i].j2i_min_send[0];
				imin[1] = xpen2real.info[i].j2i_min_send[1];

				imax[0] = xpen2real.info[i].j2i_max_send[0];
				imax[1] = xpen2real.info[i].j2i_max_send[1];

				nsend = nmap * xpen2real.info[i].j2i_ncell[0]
				             * xpen2real.info[i].j2i_ncell[1];
				nx = xpen2real.info[i].j2i_ncell_partition_send[0];
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
				for( q = imin[1]; q <= imax[1]; ++q )
				{
					for( p = imin[0]; p <= imax[0]; ++p )
					{

						if( rX )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(rhsX[pq]);
							++k;
						}
						if( rY )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(rhsY[pq]);
							++k;
						}
						if( rZ )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(rhsZ[pq]);
							++k;
						}
						if( lX )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(lhsX[pq]);
							++k;
						}
						if( lY )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(lhsY[pq]);
							++k;
						}
						if( lZ )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(lhsZ[pq]);
							++k;
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
					             * xpen2real.info[i].i2j_ncell[1];
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
					             * xpen2real.info[i].j2i_ncell[1];
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
					             * xpen2real.info[i].j2i_ncell[1];
					buffer_recv[i] = new double[ nrecv ];
					MPI_Sendrecv( buffer_send, nsend, MPI_DOUBLE, jproc, iproc, 
					              buffer_recv[i], nrecv, MPI_DOUBLE, jproc, jproc, 
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				else if( rank == jproc )
				{
					nrecv = nmap * xpen2real.info[i].i2j_ncell[0]
					             * xpen2real.info[i].i2j_ncell[1];
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
// Clear arrays
//----------------------------------------------------------------------------//
	nrecv = xpen2real.partition_recv[rank].ncell[0] 
	      * xpen2real.partition_recv[rank].ncell[1];

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

			imax[0] = xpen2real.info[i].i2j_max_recv[0];
			imax[1] = xpen2real.info[i].i2j_max_recv[1];

			nx = xpen2real.info[i].i2j_ncell_partition_recv[0];
		}
		else if( rank == iproc && xpen2real.info[i].nway > 1 )
		{
			imin[0] = xpen2real.info[i].j2i_min_recv[0];
			imin[1] = xpen2real.info[i].j2i_min_recv[1];

			imax[0] = xpen2real.info[i].j2i_max_recv[0];
			imax[1] = xpen2real.info[i].j2i_max_recv[1];

			nx = xpen2real.info[i].j2i_ncell_partition_recv[0];
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
			for( q = imin[1]; q <= imax[1]; ++q )
			{
				for( p = imin[0]; p <= imax[0]; ++p )
				{

					if( rX )
					{
						pq = q * nx + p;
						real_rhsX[pq] = buffer_recv[i][k];
						++k;
					}
					if( rY )
					{
						pq = q * nx + p;
						real_rhsY[pq] = buffer_recv[i][k];
						++k;
					}
					if( rZ )
					{
						pq = q * nx + p;
						real_rhsZ[pq] = buffer_recv[i][k];
						++k;
					}
					if( lX )
					{
						pq = q * nx + p;
						real_lhsX[pq] = buffer_recv[i][k];
						++k;
					}
					if( lY )
					{
						pq = q * nx + p;
						real_lhsY[pq] = buffer_recv[i][k];
						++k;
					}
					if( lZ )
					{
						pq = q * nx + p;
						real_lhsZ[pq] = buffer_recv[i][k];
						++k;
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
