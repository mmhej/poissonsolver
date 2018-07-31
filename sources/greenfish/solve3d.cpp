//----------------------------------------------------------------------------//
/*
  File:         solve3d.cpp

  Description:  
*/
//----------------------------------------------------------------------------//
void class_greenfish::solve3d(  )
{

//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//
	const double pi = acos(-1.0);

//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int     nproc, rank;
	int     i, j, k, kn, kjn, ijk;
	int     nfft;
	int     ncell[3];

	std::complex<double> div_psi;

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
			lZ = true;
		}
		else
		{
			std::cerr << " [greenfish.solve2d]: "
			          << "Error rhs fields have not been pushed correctly for grad option" 
			          << std::endl;
			return;
		}
	}
	else if(lhs_div)
	{
		if(rX && rY && rZ)
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
	else if(lhs_curl)
	{
		if(rX && rY && rZ)
		{
			lX = true;
			lY = true;
			lZ = true;
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
	ncell[2] = xpen[rank].ncell[2];

	pen_rhs.resize( ncell[0] );

	for (k = 0; k < ncell[2]; ++k )
	{
		kn = k * ncell[1];
		for (j = 0; j < ncell[1]; ++j )
		{
			kjn = (kn + j) * ncell[0];
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store mesh array in x-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for(i = 0; i < ncell[0]; ++i )
			{
				ijk = kjn + i;
				if( rX )
				{
					pen_rhs.X[i] = rhsX[ijk];
				}
				if( rY )
				{
					pen_rhs.Y[i] = rhsY[ijk];
				}
				if( rZ )
				{
					pen_rhs.Z[i] = rhsZ[ijk];
				}
			}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Fourier transform pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			pen_rhs.fft( );
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store pencil in fft array
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for (i = 0; i < ncell[0]; ++i )
			{
				ijk = kjn + i;
				if( rX )
				{
					rhsX[ijk] = pen_rhs.X[i];
				}
				if( rY )
				{
					rhsY[ijk] = pen_rhs.Y[i];
				}
				if( rZ )
				{
					rhsZ[ijk] = pen_rhs.Z[i];
				}
			}
		}
	}

//----------------------------------------------------------------------------//
// Map to y-pencil partition
//----------------------------------------------------------------------------//
	map( xpen2ypen );

//----------------------------------------------------------------------------//
// FFT y-pencils
//----------------------------------------------------------------------------//
	ncell[0] = ypen[rank].ncell[0];
	ncell[1] = ypen[rank].ncell[1];
	ncell[2] = ypen[rank].ncell[2];

	pen_rhs.resize( ncell[1] );

	for (k = 0; k < ncell[2]; ++k )
	{
		kn = k * ncell[1];
		for (i = 0; i < ncell[0]; ++i )
		{
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store mesh array in y-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for(j = 0; j < ncell[1]; ++j )
			{
				ijk = (kn + j) * ncell[0] + i;
				if( rX )
				{
					pen_rhs.X[j] = rhsX[ijk];
				}
				if( rY )
				{
					pen_rhs.Y[j] = rhsY[ijk];
				}
				if( rZ )
				{
					pen_rhs.Z[j] = rhsZ[ijk];
				}
			}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Fourier transform pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			pen_rhs.fft( );
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store pencil in fft array
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for (j = 0; j < ncell[1]; ++j )
			{
				ijk = (kn + j) * ncell[0] + i;
				if( rX )
				{
					rhsX[ijk] = pen_rhs.X[j];
				}
				if( rY )
				{
					rhsY[ijk] = pen_rhs.Y[j];
				}
				if( rZ )
				{
					rhsZ[ijk] = pen_rhs.Z[j];
				}
			}
		}
	}


//----------------------------------------------------------------------------//
// FFT z-pencils and perform Fourier space operations
//----------------------------------------------------------------------------//
	if(domain_bc[2] == 0)
	{
		nfft = 2*ncell[2];
	}
	else
	{
		nfft = ncell[2];
	}

	pen_rhs.resize( nfft );
	pen_lhs.resize( nfft );

	if( lX )
	{
		lhsX = new std::complex<double>[ncell[0]*ncell[1]*ncell[2]]();
	}
	if( lY )
	{
		lhsY = new std::complex<double>[ncell[0]*ncell[1]*ncell[2]]();
	}
	if( lZ )
	{
		lhsZ = new std::complex<double>[ncell[0]*ncell[1]*ncell[2]]();
	}


	for (i = 0; i < ncell[0]; ++i )
	{
		for(j = 0; j < ncell[1]; ++j )
		{

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store mesh array in z-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for(k = 0; k < ncell[2]; ++k )
			{
				ijk = (k * ncell[1] + j) * ncell[0] + i;

				if( rX )
				{
					pen_rhs.X[k] = rhsX[ijk];
				}
				if( rY )
				{
					pen_rhs.Y[k] = rhsY[ijk];
				}
				if( rZ )
				{
					pen_rhs.Z[k] = rhsZ[ijk];
				}
			}
// Zero-padding
			if( nfft > ncell[2] )
			{
				for(k = ncell[2]; k < nfft; ++k )
				{
					if( rX )
					{
						pen_rhs.X[k] = {0.0,0.0};
					}
					if( rY )
					{
						pen_rhs.Y[k] = {0.0,0.0};
					}
					if( rZ )
					{
						pen_rhs.Z[k] = {0.0,0.0};
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
			for (k = 0; k < nfft; ++k )
			{
				ijk = (k * ncell[1] + j) * ncell[0] + i;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Regularise rhs
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				if( regularisation > 0 )
				{
					if( rX )
					{
						pen_rhs.X[k] *= zeta[ijk];
					}
					if( rY )
					{
						pen_rhs.Y[k] *= zeta[ijk];
					}
					if( rZ )
					{
						pen_rhs.Z[k] *= zeta[ijk];
					}
				}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Do convolution
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				if(lhs_grad)
				{
					pen_lhs.X[k] = ikX[i] * G3D[ijk] * pen_rhs.X[k];
					pen_lhs.Y[k] = ikY[j] * G3D[ijk] * pen_rhs.Y[k];
					pen_lhs.Z[k] = ikZ[k] * G3D[ijk] * pen_rhs.Z[k];
				}
				else if(lhs_div)
				{
					pen_lhs.X[k] = G3D[ijk] * ( ikX[i] * pen_rhs.X[k]
					                          + ikY[j] * pen_rhs.Y[k]
					                          + ikZ[k] * pen_rhs.Z[k] );
				}
				else if(lhs_curl)
				{
					pen_lhs.X[k] = G3D[ijk] * ( ikY[j] * pen_rhs.Z[k]
					                          - ikZ[k] * pen_rhs.Y[k] );
					pen_lhs.Y[k] = G3D[ijk] * ( ikZ[k] * pen_rhs.X[k]
					                          - ikX[i] * pen_rhs.Z[k] );
					pen_lhs.Z[k] = G3D[ijk] * ( ikX[i] * pen_rhs.Y[k]
					                          - ikY[j] * pen_rhs.X[k] );
				}
				else
				{
					if( lX )
					{
//						pen_lhs.X[k] = pen_rhs.X[k];
						pen_lhs.X[k] = G3D[ijk] * pen_rhs.X[k];
					}
					if( lY )
					{
//						pen_lhs.Y[k] = pen_rhs.Y[k];
						pen_lhs.Y[k] = G3D[ijk] * pen_rhs.Y[k];
					}
					if( lZ )
					{
//						pen_lhs.Z[k] = pen_rhs.Z[k];
						pen_lhs.Z[k] = G3D[ijk] * pen_rhs.Z[k];
					}
				}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Reproject rhs field
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				if( rhs_reproject && rX && rY && rZ )
				{
					div_psi = G3D[ijk] * ( ikX[i] * pen_rhs.X[k]
					                     + ikY[j] * pen_rhs.Y[k]
					                     + ikZ[k] * pen_rhs.Z[k] );
					pen_rhs.X[k] += ikX[i] * div_psi;
					pen_rhs.Y[k] += ikY[j] * div_psi;
					pen_rhs.Z[k] += ikZ[k] * div_psi;
				}

			}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// IFFT z-pencils
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			pen_rhs.ifft( );
			pen_lhs.ifft( );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store solution
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for (k = 0; k < ncell[2]; ++k )
			{
				ijk = (k * ncell[1] + j) * ncell[0] + i;

				if( rX )
				{
					rhsX[ijk] = pen_rhs.X[k]/(double)nfft;
				}
				if( rY )
				{
					rhsY[ijk] = pen_rhs.Y[k]/(double)nfft;
				}
				if( rZ )
				{
					rhsZ[ijk] = pen_rhs.Z[k]/(double)nfft;
				}

				if( lX )
				{
					lhsX[ijk] = pen_lhs.X[k]/(double)nfft;
				}
				if( lY )
				{
					lhsY[ijk] = pen_lhs.Y[k]/(double)nfft;
				}
				if( lZ )
				{
					lhsZ[ijk] = pen_lhs.Z[k]/(double)nfft;
				}

			}
		}
	}


//----------------------------------------------------------------------------//
// IFFT y-pencils
//----------------------------------------------------------------------------//
	pen_rhs.resize( ncell[1] );
	pen_lhs.resize( ncell[1] );

	for (k = 0; k < ncell[2]; ++k )
	{
		kn = k * ncell[1];
		for (i = 0; i < ncell[0]; ++i )
		{
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store mesh array in y-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for(j = 0; j < ncell[1]; ++j )
			{
				ijk = (kn + j) * ncell[0] + i;

				if( rX )
				{
					pen_rhs.X[j] = rhsX[ijk];
				}
				if( rY )
				{
					pen_rhs.Y[j] = rhsY[ijk];
				}
				if( rZ )
				{
					pen_rhs.Z[j] = rhsZ[ijk];
				}

				if( lX )
				{
					pen_lhs.X[j] = lhsX[ijk];
				}
				if( lY )
				{
					pen_lhs.Y[j] = lhsY[ijk];
				}
				if( lZ )
				{
					pen_lhs.Z[j] = lhsZ[ijk];
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
			for(j = 0; j < ncell[1]; ++j )
			{
				ijk = (kn + j) * ncell[0] + i;

				if( rX )
				{
					rhsX[ijk] = pen_rhs.X[j]/(double)ncell[1];
				}
				if( rY )
				{
					rhsY[ijk] = pen_rhs.Y[j]/(double)ncell[1];
				}
				if( rZ )
				{
					rhsZ[ijk] = pen_rhs.Z[j]/(double)ncell[1];
				}

				if( lX )
				{
					lhsX[ijk] = pen_lhs.X[j]/(double)ncell[1];
				}
				if( lY )
				{
					lhsY[ijk] = pen_lhs.Y[j]/(double)ncell[1];
				}
				if( lZ )
				{
					lhsZ[ijk] = pen_lhs.Z[j]/(double)ncell[1];
				}
			}

		}
	}


//----------------------------------------------------------------------------//
// Map to x-pencils
//----------------------------------------------------------------------------//
	map( ypen2xpen );

//----------------------------------------------------------------------------//
// IFFT x-pencils
//----------------------------------------------------------------------------//
	ncell[0] = xpen[rank].ncell[0];
	ncell[1] = xpen[rank].ncell[1];
	ncell[2] = xpen[rank].ncell[2];

	pen_rhs.resize( ncell[0] );
	pen_lhs.resize( ncell[0] );

	for (k = 0; k < ncell[2]; ++k )
	{
		kn = k * ncell[1];
		for (j = 0; j < ncell[1]; ++j )
		{
			kjn = (kn + j) * ncell[0];
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store mesh array in x-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for(i = 0; i < ncell[0]; ++i )
			{
				ijk = kjn + i;

				if( rX )
				{
					pen_rhs.X[i] = rhsX[ijk];
				}
				if( rY )
				{
					pen_rhs.Y[i] = rhsY[ijk];
				}
				if( rZ )
				{
					pen_rhs.Z[i] = rhsZ[ijk];
				}

				if( lX )
				{
					pen_lhs.X[i] = lhsX[ijk];
				}
				if( lY )
				{
					pen_lhs.Y[i] = lhsY[ijk];
				}
				if( lZ )
				{
					pen_lhs.Z[i] = lhsZ[ijk];
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
				ijk = kjn + i;

				if( rX )
				{
					rhsX[ijk] = pen_rhs.X[i]/(double)ncell[0];
				}
				if( rY )
				{
					rhsY[ijk] = pen_rhs.Y[i]/(double)ncell[0];
				}
				if( rZ )
				{
					rhsZ[ijk] = pen_rhs.Z[i]/(double)ncell[0];
				}

				if( lX )
				{
					lhsX[ijk] = pen_lhs.X[i]/(double)ncell[0];
				}
				if( lY )
				{
					lhsY[ijk] = pen_lhs.Y[i]/(double)ncell[0];
				}
				if( lZ )
				{
					lhsZ[ijk] = pen_lhs.Z[i]/(double)ncell[0];
				}
			}

		}
	}

//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}


