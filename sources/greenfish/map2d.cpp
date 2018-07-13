//----------------------------------------------------------------------------//
/*
  File:         map2d.cpp

  Description:  
*/
//----------------------------------------------------------------------------//

void class_greenfish::map2d( class_communication comm )
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

	bool mG = false;

//----------------------------------------------------------------------------//
// Objects
//----------------------------------------------------------------------------//
	double * buffer_send = NULL;
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

//		if(rank == 0){ std::cout << "COMM: " << i << std::endl; }
//		MPI_Barrier(MPI_COMM_WORLD);
//----------------------------------------------------------------------------//
// Self-communication
//----------------------------------------------------------------------------//
		if( comm.info[i].nway == 0 )
		{

			imin[0] = comm.info[i].i2j_min_send[0];
			imin[1] = comm.info[i].i2j_min_send[1];

			imax[0] = comm.info[i].i2j_max_send[0];
			imax[1] = comm.info[i].i2j_max_send[1];

			nx = comm.info[i].i2j_ncell_partition_send[0];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Fill receive buffer directly
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			nsend = 2 * nmap * comm.info[i].i2j_ncell[0]
			                 * comm.info[i].i2j_ncell[1];
			buffer_recv[i] = new double[ nsend ]();

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
						buffer_recv[i][k] = std::imag(rhsX[pq]);
						++k;
					}
					if( rY )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(rhsY[pq]);
						++k;
						buffer_recv[i][k] = std::imag(rhsY[pq]);
						++k;
					}
					if( rZ )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(rhsZ[pq]);
						++k;
						buffer_recv[i][k] = std::imag(rhsZ[pq]);
						++k;
					}

					if( lX )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(lhsX[pq]);
						++k;
						buffer_recv[i][k] = std::imag(lhsX[pq]);
						++k;
					}
					if( lY )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(lhsY[pq]);
						++k;
						buffer_recv[i][k] = std::imag(lhsY[pq]);
						++k;
					}
					if( lZ )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(lhsZ[pq]);
						++k;
						buffer_recv[i][k] = std::imag(lhsZ[pq]);
						++k;
					}

					if( mG )
					{
						pq = q * nx + p;
						buffer_recv[i][k] = std::real(mapG[pq]);
						++k;
						buffer_recv[i][k] = std::imag(mapG[pq]);
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

				imax[0] = comm.info[i].i2j_max_send[0];
				imax[1] = comm.info[i].i2j_max_send[1];

				nsend = 2 * nmap * comm.info[i].i2j_ncell[0]
				                 * comm.info[i].i2j_ncell[1];
				nx = comm.info[i].i2j_ncell_partition_send[0];
			}
			else if( rank == jproc && comm.info[i].nway > 1 )
			{

				imin[0] = comm.info[i].j2i_min_send[0];
				imin[1] = comm.info[i].j2i_min_send[1];

				imax[0] = comm.info[i].j2i_max_send[0];
				imax[1] = comm.info[i].j2i_max_send[1];

				nsend = 2 * nmap * comm.info[i].j2i_ncell[0]
				                 * comm.info[i].j2i_ncell[1];
				nx = comm.info[i].j2i_ncell_partition_send[0];
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
				for( q = imin[1]; q <= imax[1]; ++q )
				{
					for( p = imin[0]; p <= imax[0]; ++p )
					{

						if( rX )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(rhsX[pq]);
							++k;
							buffer_send[k] = std::imag(rhsX[pq]);
							++k;
						}
						if( rY )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(rhsY[pq]);
							++k;
							buffer_send[k] = std::imag(rhsY[pq]);
							++k;
						}
						if( rZ )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(rhsZ[pq]);
							++k;
							buffer_send[k] = std::imag(rhsZ[pq]);
							++k;
						}
						if( lX )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(lhsX[pq]);
							++k;
							buffer_send[k] = std::imag(lhsX[pq]);
							++k;
						}
						if( lY )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(lhsY[pq]);
							++k;
							buffer_send[k] = std::imag(lhsY[pq]);
							++k;
						}
						if( lZ )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(lhsZ[pq]);
							++k;
							buffer_send[k] = std::imag(lhsZ[pq]);
							++k;
						}

						if( mG )
						{
							pq = q * nx + p;
							buffer_send[k] = std::real(mapG[pq]);
							++k;
							buffer_send[k] = std::imag(mapG[pq]);
							++k;
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
					MPI_Send( buffer_send, nsend, MPI_C_DOUBLE_COMPLEX, jproc, jproc,
					          MPI_COMM_WORLD );
				}
				else if( rank == jproc ) // Receiver
				{
					nrecv = 2 * nmap * comm.info[i].i2j_ncell[0] 
					                 * comm.info[i].i2j_ncell[1];
					buffer_recv[i] = new double[ nrecv ]();
					MPI_Recv( buffer_recv[i], nrecv, MPI_DOUBLE, iproc, jproc, 
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
			}
			else if( comm.info[i].nway == 2 ) // only j to i
			{
				if( rank == jproc ) // Sender
				{
					MPI_Send( buffer_send, nsend, MPI_C_DOUBLE_COMPLEX, iproc, iproc,
					          MPI_COMM_WORLD );
				}
				else if( rank == iproc ) // Receiver
				{
					nrecv = 2 * nmap * comm.info[i].j2i_ncell[0] 
					                 * comm.info[i].j2i_ncell[1];
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
					                 * comm.info[i].j2i_ncell[1];
					buffer_recv[i] = new double[ nrecv ]();
					MPI_Sendrecv( buffer_send, nsend, MPI_DOUBLE, jproc, iproc, 
					              buffer_recv[i], nrecv, MPI_DOUBLE, jproc, jproc, 
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				else if( rank == jproc )
				{
					nrecv = 2 * nmap * comm.info[i].i2j_ncell[0]
					                 * comm.info[i].i2j_ncell[1];
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
	      * comm.partition_recv[rank].ncell[1];
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

			imax[0] = comm.info[i].i2j_max_recv[0];
			imax[1] = comm.info[i].i2j_max_recv[1];

			nx = comm.info[i].i2j_ncell_partition_recv[0];
		}
		else if( rank == iproc && comm.info[i].nway > 1 )
		{
			imin[0] = comm.info[i].j2i_min_recv[0];
			imin[1] = comm.info[i].j2i_min_recv[1];

			imax[0] = comm.info[i].j2i_max_recv[0];
			imax[1] = comm.info[i].j2i_max_recv[1];

			nx = comm.info[i].j2i_ncell_partition_recv[0];
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
						rhsX[pq] = { buffer_recv[i][k], buffer_recv[i][k+1] };
						k += 2;
					}
					if( rY )
					{
						pq = q * nx + p;
						rhsY[pq] = { buffer_recv[i][k], buffer_recv[i][k+1] };
						k += 2;
					}
					if( rZ )
					{
						pq = q * nx + p;
						rhsZ[pq] = { buffer_recv[i][k], buffer_recv[i][k+1] };
						k += 2;
					}
					if( lX )
					{
						pq = q * nx + p;
						lhsX[pq] = { buffer_recv[i][k], buffer_recv[i][k+1] };
						k += 2;
					}
					if( lY )
					{
						pq = q * nx + p;
						lhsY[pq] = { buffer_recv[i][k], buffer_recv[i][k+1] };
						k += 2;
					}
					if( lZ )
					{
						pq = q * nx + p;
						lhsZ[pq] = { buffer_recv[i][k], buffer_recv[i][k+1] };
						k += 2;
					}

					if( mG )
					{
						pq = q * nx + p;
						mapG[pq] = { buffer_recv[i][k], buffer_recv[i][k+1] };
						k += 2;
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
