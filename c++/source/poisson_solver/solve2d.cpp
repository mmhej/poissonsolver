//----------------------------------------------------------------------------//
/*
  File:         solve2d.cpp

  Description:  Solves the Poisson equation in 2D
*/
//----------------------------------------------------------------------------//
void poisson_solver::solve2d(  )
{

//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//
	const double pi = 3.1415926535897932;

//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int     nproc, rank;
	int     i, j, jn, ij;
	int     nfft;
	int     ncell[2];

	bool rX = false;
	bool rY = false;
	bool rZ = false;
	bool lX = false;
	bool lY = false;
	bool lZ = false;

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//----------------------------------------------------------------------------//
// Find out which fields to map
//----------------------------------------------------------------------------//
	if( rhsX != NULL ){ rX = true; }
	if( rhsY != NULL ){ rY = true; }
	if( rhsZ != NULL ){ rZ = true; }

	if(lhs_grad)
	{
		if(rX)
		{
			lX = true;
			lY = true;
		}
		else
		{
			std::cerr << " [poisson_solver.solve2d]: "
			          << "Error rhs fields have not been pushed correctly for grad option" 
			          << std::endl;
			return;
		}
	}
	else if(lhs_div)
	{
		if(rX && rY)
		{
			lX = true;
		}
		else
		{
			std::cerr << " [poisson_solver.solve2d]: "
			          << "Error rhs fields have not been pushed correctly for div option" 
			          << std::endl;
			return;
		}
	}
	else if(lhs_curl)
	{
		if(rZ)
		{
			lX = true;
			lY = true;
		}
		else
		{
			std::cerr << " [poisson_solver.solve2d]: "
			          << "Error rhs fields have not been pushed correctly for curl option" 
			          << std::endl;
			return;
		}
	}
	else
	{
		if(rX)
		{
			lX = true;
		}
		if(rY)
		{
			lY = true;
		}
		if(rZ)
		{
			lZ = true;
		}
	}

//----------------------------------------------------------------------------//
// Construct pencils
//----------------------------------------------------------------------------//
	class_pencil pen_rhs(rX,rY,rZ);
	class_pencil pen_lhs(lX,lY,lZ);

//----------------------------------------------------------------------------//
// FFT x-pencils
//----------------------------------------------------------------------------//
	ncell[0] = xpen[rank].ncell[0];
	ncell[1] = xpen[rank].ncell[1];

	pen_rhs.resize( ncell[0] );

	for (j = 0; j < ncell[1]; ++j )
	{
		jn = j * ncell[0];
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft mesh array in x-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		for(i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			if( rX )
			{
				pen_rhs.X[i] = rhsX[ij];
			}
			if( rY )
			{
				pen_rhs.Y[i] = rhsY[ij];
			}
			if( rZ )
			{
				pen_rhs.Z[i] = rhsZ[ij];
			}
		}
		pen_rhs.fft( );
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			if( rX )
			{
				rhsX[ij] = pen_rhs.X[i];
			}
			if( rY )
			{
				rhsY[ij] = pen_rhs.Y[i];
			}
			if( rZ )
			{
				rhsZ[ij] = pen_rhs.Z[i];
			}
		}
	}

//----------------------------------------------------------------------------//
// Map to y-pencils
//----------------------------------------------------------------------------//
	map( xpen2ypen );

/*
	if(rank == 0)
	{
		for(i = 0; i < 3*16; i++ )
		{
			std::cout << i << " " << rhsX[i] << std::endl;
		}
	}
*/

//----------------------------------------------------------------------------//
// FFT y-pencils and perform Fourier space operations
//----------------------------------------------------------------------------//
	ncell[0] = ypen[rank].ncell[0];
	ncell[1] = ypen[rank].ncell[1];

	if(domain_bc[1] == 0)
	{
		nfft = 2*ncell[1];
	}
	else
	{
		nfft = ncell[1];
	}

	pen_rhs.resize( nfft );
	pen_lhs.resize( nfft );

	if( lX )
	{
		lhsX = new std::complex<double>[ncell[0]*ncell[1]]();
	}
	if( lY )
	{
		lhsY = new std::complex<double>[ncell[0]*ncell[1]]();
	}
	if( lZ )
	{
		lhsZ = new std::complex<double>[ncell[0]*ncell[1]]();
	}


	for (i = 0; i < ncell[0]; ++i )
	{
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft mesh array in y-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		for(j = 0; j < ncell[1]; ++j )
		{
			ij = j * ncell[0] + i;
			if( rX )
			{
				pen_rhs.X[j] = rhsX[ij];
			}
			if( rY )
			{
				pen_rhs.Y[j] = rhsY[ij];
			}
			if( rZ )
			{
				pen_rhs.Z[j] = rhsZ[ij];
			}
		}
// Zero-padding
		if( nfft > ncell[1] )
		{
			for(j = ncell[1]; j < nfft; ++j )
			{
				if( rX )
				{
					pen_rhs.X[j] = {0.0,0.0};
				}
				if( rY )
				{
					pen_rhs.Y[j] = {0.0,0.0};
				}
				if( rZ )
				{
					pen_rhs.Z[j] = {0.0,0.0};
				}
			}
		}

/*
if(rank == 0 && i == 0)
{
	for(j = 0; j < ncell[1]; j++ )
	{
		ij = j * ncell[0] + i;
		std::cout << rhsX[j] << " " << rhsX[ij] << std::endl;
	}
}
*/

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// FFT y-pencils
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		pen_rhs.fft( );

//----------------------------------------------------------------------------//
// Perform Fourier space operations
//----------------------------------------------------------------------------//
		for (j = 0; j < nfft; ++j )
		{
			ij = j * ncell[0] + i;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Regularise rhs
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			if( regularisation > 0 )
			{
				if( rX )
				{
					pen_rhs.X[j] *= zeta[ij];
				}
				if( rY )
				{
					pen_rhs.Y[j] *= zeta[ij];
				}
				if( rZ )
				{
					pen_rhs.Z[j] *= zeta[ij];
				}
			}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Do convolution
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			if(lhs_grad)
			{

//				pen_lhs.X[j] = ikX[i] * pen_rhs.X[j];
//				pen_lhs.Y[j] = ikY[j] * pen_rhs.X[j];

				pen_lhs.X[j] = ikX[i] * G2D[ij] * pen_rhs.X[j];
				pen_lhs.Y[j] = ikY[j] * G2D[ij] * pen_rhs.X[j];
			}
			else if(lhs_div)
			{
				pen_lhs.X[j] = G2D[ij] * ( ikX[i] * pen_rhs.X[j]
				                         + ikY[j] * pen_rhs.Y[j] );
			}
			else if(lhs_curl)
			{
				pen_lhs.X[j] =   ikY[j] * G2D[ij] * pen_rhs.Z[j];
				pen_lhs.Y[j] = - ikX[i] * G2D[ij] * pen_rhs.Z[j];
			}
			else
			{
				if( lX )
				{
					pen_lhs.X[j] = G2D[ij] * pen_rhs.X[j];
				}
				if( lY )
				{
					pen_lhs.Y[j] = G2D[ij] * pen_rhs.Y[j];
				}
				if( lZ )
				{
					pen_lhs.Z[j] = G2D[ij] * pen_rhs.Z[j];
				}
			}
		}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// IFFT y-pencils
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		pen_rhs.ifft( );
		pen_lhs.ifft( );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store solution
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		for (j = 0; j < ncell[1]; ++j )
		{
			ij = j * ncell[0] + i;
			if( rX )
			{
				rhsX[ij] = pen_rhs.X[j]/(double)nfft;
			}
			if( rY )
			{
				rhsY[ij] = pen_rhs.Y[j]/(double)nfft;
			}
			if( rZ )
			{
				rhsZ[ij] = pen_rhs.Z[j]/(double)nfft;
			}

			if( lX )
			{
				lhsX[ij] = pen_lhs.X[j]/(double)nfft;
			}
			if( lY )
			{
				lhsY[ij] = pen_lhs.Y[j]/(double)nfft;
			}
			if( lZ )
			{
				lhsZ[ij] = pen_lhs.Z[j]/(double)nfft;
			}
		}
	}


//	if(rank == 0)
//	{
//		for(i = 0; i < 10; i++ )
//		{
//			std::cout << lhsX[i] << std::endl;
//		}
//	}


//----------------------------------------------------------------------------//
// Map to x-pencils
//----------------------------------------------------------------------------//
	map( ypen2xpen );

//----------------------------------------------------------------------------//
// IFFT x-pencils
//----------------------------------------------------------------------------//
	ncell[0] = xpen[rank].ncell[0];
	ncell[1] = xpen[rank].ncell[1];

	pen_rhs.resize( ncell[0] );
	pen_lhs.resize( ncell[0] );

	for (j = 0; j < ncell[1]; ++j )
	{
		jn = j * ncell[0];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft mesh array in x-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			if( rX )
			{
				pen_rhs.X[i] = rhsX[ij];
			}
			if( rY )
			{
				pen_rhs.Y[i] = rhsY[ij];
			}
			if( rZ )
			{
				pen_rhs.Z[i] = rhsZ[ij];
			}

			if( lX )
			{
				pen_lhs.X[i] = lhsX[ij];
			}
			if( lY )
			{
				pen_lhs.Y[i] = lhsY[ij];
			}
			if( lZ )
			{
				pen_lhs.Z[i] = lhsZ[ij];
			}
		}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Perform ifft
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		pen_rhs.ifft( );
		pen_lhs.ifft( );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft field
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			if( rX )
			{
				rhsX[ij] = pen_rhs.X[i]/(double)ncell[0];
			}
			if( rY )
			{
				rhsY[ij] = pen_rhs.Y[i]/(double)ncell[0];
			}
			if( rZ )
			{
				rhsZ[ij] = pen_rhs.Z[i]/(double)ncell[0];
			}

			if( lX )
			{
				lhsX[ij] = pen_lhs.X[i]/(double)ncell[0];
			}
			if( lY )
			{
				lhsY[ij] = pen_lhs.Y[i]/(double)ncell[0];
			}
			if( lZ )
			{
				lhsZ[ij] = pen_lhs.Z[i]/(double)ncell[0];
			}
		}
	}

//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}


