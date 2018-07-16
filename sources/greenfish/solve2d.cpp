//----------------------------------------------------------------------------//
/*
  File:         solve2d.cpp

  Description:  
*/
//----------------------------------------------------------------------------//
void class_greenfish::solve2d(  )
{

//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//
	const double pi = acos(-1.0);
	const std::complex<double> im = {0.0,1.0};

//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int     nproc, rank;
	int     i, j, ij;
	int     nfft;
	int     ncell[2];
	int     icell[2];
	double  dx[2];

	std::complex<double> G;

	double dk;
	double * kX;
	double * kY;

	double kL;
	double kx2, ky2;

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

	if(rhs_grad)
	{
std::cout << "Solving grad" << std::endl;
		if(rX)
		{
			lX = true;
			lY = true;
		}
		else
		{
			std::cerr << " [greenfish.solve2d]: "
			          << "Error rhs fields have not been pushed correctly for grad option" 
			          << std::endl;
			return;
		}
	}
	if(rhs_div)
	{
std::cout << "Solving div" << std::endl;
		if(rX && rY)
		{
			lX = true;
		}
		else
		{
			std::cerr << " [greenfish.solve2d]: "
			          << "Error rhs fields have not been pushed correctly for div option" 
			          << std::endl;
			return;
		}
	}
	else if(rhs_curl)
	{
std::cout << "Solving curl" << std::endl;
		if(rZ)
		{
			lX = true;
			lY = true;
		}
		else
		{
			std::cerr << " [greenfish.solve2d]: "
			          << "Error rhs fields have not been pushed correctly for curl option" 
			          << std::endl;
			return;
		}
	}
	else
	{
std::cout << "Solving regular" << std::endl;
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
	ncell[0] = partition.xpen[rank].ncell[0];
	ncell[1] = partition.xpen[rank].ncell[1];

	pen_rhs.resize( ncell[0] );

	for (j = 0; j < ncell[1]; ++j )
	{

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft mesh array in x-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		for(i = 0; i < ncell[0]; ++i )
		{
			ij = j * ncell[0] + i;
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
			ij = j * ncell[0] + i;
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

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Calculate wavenumbers in x-dir and store in global array
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	dk = 1.0/double(ncell[0]);
	kX = new double[ncell[0]]();

	for (i = 0; i < ncell[0]; ++i )
	{

		if( 2*i < ncell[0] )
		{
			kX[i] = double(i) * dk;
		}
		else
		{
			kX[i] = double(i) * dk - 1.0;
		}
//		if(rank == 0){ std::cout << kX[i] << std::endl; }
	}

//----------------------------------------------------------------------------//
// Map to y-pencils
//----------------------------------------------------------------------------//
	map2d( xpen2ypen );

//----------------------------------------------------------------------------//
// FFT y-pencils and perform Fourier space operations
//----------------------------------------------------------------------------//
	ncell[0] = partition.ypen[rank].ncell[0];
	ncell[1] = partition.ypen[rank].ncell[1];

	icell[0] = partition.ypen[rank].icell[0];
	icell[1] = partition.ypen[rank].icell[1];

	dx[0]    = partition.ypen[rank].dx[0];
	dx[1]    = partition.ypen[rank].dx[1];

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

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Calculate wavenumbers in y-dir and store in global array
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
//	dk = 1.0/double(ncell[1]);
//	kY = new double[ncell[1]];
/*
	dk = 1.0/double(nfft);
	kY = new double[nfft]();

//	for (j = 0; j < ncell[1]; ++j )
	for (j = 0; j < nfft; ++j )
	{
//		if( 2*j < ncell[1] )
		if( 2*j < nfft )
		{
			kY[j] = double(j) * dk;
		}
		else
		{
			kY[j] = double(j) * dk - 1.0;
		}
	}
*/

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

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// FFT y-pencils
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		pen_rhs.fft( );

//----------------------------------------------------------------------------//
// Perform Fourier space operations
//----------------------------------------------------------------------------//
//		for (j = 0; j < ncell[1]; ++j )
		for (j = 0; j < nfft; ++j )
		{
			ij = j * ncell[0] + i;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Do convolution
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			if(rhs_grad)
			{
				pen_lhs.X[j] = - ikX[i] * rhsG[ij] * pen_rhs.X[j];
				pen_lhs.Y[j] = - ikY[j] * rhsG[ij] * pen_rhs.X[j];
			}
			if(rhs_div)
			{
				pen_lhs.X[j] = - ikX[i] * rhsG[ij] * pen_rhs.X[j];
				               - ikY[j] * rhsG[ij] * pen_rhs.Y[j];
			}
			else if(rhs_curl)
			{
				pen_lhs.X[j] =   ikY[j] * rhsG[ij] * pen_rhs.Z[j];
				pen_lhs.Y[j] = - ikX[i] * rhsG[ij] * pen_rhs.Z[j];
			}
			else
			{
				if( lX )
				{
					pen_lhs.X[j] = rhsG[ij] * pen_rhs.X[j];
				}
				if( lY )
				{
					pen_lhs.Y[j] = rhsG[ij] * pen_rhs.Y[j];
				}
				if( lZ )
				{
					pen_lhs.Z[j] = rhsG[ij] * pen_rhs.Z[j];
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

//----------------------------------------------------------------------------//
// Map to x-pencils
//----------------------------------------------------------------------------//
	map2d( ypen2xpen );

//----------------------------------------------------------------------------//
// IFFT x-pencils
//----------------------------------------------------------------------------//
	ncell[0] = partition.xpen[rank].ncell[0];
	ncell[1] = partition.xpen[rank].ncell[1];

	pen_rhs.resize( ncell[0] );
	pen_lhs.resize( ncell[0] );

	for (j = 0; j < ncell[1]; ++j )
	{

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft mesh array in x-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = j * ncell[0] + i;
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
			ij = j * ncell[0] + i;
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


