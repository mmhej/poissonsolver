//----------------------------------------------------------------------------//
/*
  File:         greens3d.cpp

  Description:  
*/
//----------------------------------------------------------------------------//
void class_greenfish::greens3d(  )
{

//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//
	const double pi = acos(-1.0);

	const double c_1_2pi2 = 0.5/pow(pi,2);

//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int     nproc, rank;
	int     i, j, k, kn, kjn, ijk;
	int     nfft;
	int     ncell[3];
	int     icell[3];
	double  dx[3];
	double  xmin[3];

	double  x,y,z,r,rho;
	double  sigma;
	double  dk;

	double  C;

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//============================================================================//
// Create extended partitions and communications
//============================================================================//
//	class_communication comm;
	std::vector<class_partition> xpen_ext;
	std::vector<class_partition> ypen_ext;

	dx[0] = domain_dx[0];
	dx[1] = domain_dx[1];
	dx[2] = domain_dx[2];

	xpen_ext = partition_setup( 0, domain_ncell, domain_bc, domain_dx, true, true, true );
	ypen_ext = partition_setup( 1, domain_ncell, domain_bc, domain_dx, true, true, true );

	class_communication comm( xpen_ext, ypen_ext );

//============================================================================//
// Create spectral derivative operator on y-pencil partition
//============================================================================//
	ncell[0] = ypen_ext[rank].ncell[0];
	ncell[1] = ypen_ext[rank].ncell[1];
	ncell[2] = ypen_ext[rank].ncell[2];

	icell[0] = ypen_ext[rank].icell[0];
	icell[1] = ypen_ext[rank].icell[1];
	icell[2] = ypen_ext[rank].icell[2];

	nfft = xpen_ext[rank].ncell[0];

	dk  = 1.0/( double( nfft ) * dx[0] );
	ikX = new std::complex<double>[ncell[0]]();

	for(i = 0; i < ncell[0]; ++i )
	{
		if( 2*(icell[0]+i) < nfft )
		{
			ikX[i] = {0.0, 2.0 * pi * double(icell[0]+i) * dk};
		}
		else
		{
			ikX[i] = {0.0, 2.0 * pi * double(icell[0]+i-nfft) * dk };
		}
	}

	dk  = 1.0/( double(ncell[1]) * dx[1] );
	ikY = new std::complex<double>[ncell[1]]();

	for(j = 0; j < ncell[1]; ++j )
	{
		if( 2*j < ncell[1] )
		{
			ikY[j] = {0.0, 2.0 * pi * double(j) * dk};
		}
		else
		{
			ikY[j] = {0.0, 2.0 * pi * double(j-ncell[1]) * dk};
		}
	}

	dk  = 1.0/( double(ncell[2]) * dx[2] );
	ikZ = new std::complex<double>[ncell[2]]();

	for(k = 0; k < ncell[2]; ++k )
	{
		if( 2*k < ncell[2] )
		{
			ikZ[k] = {0.0, 2.0 * pi * double(k) * dk};
		}
		else
		{
			ikZ[k] = {0.0, 2.0 * pi * double(k-ncell[2]) * dk};
		}
	}

//============================================================================//
// Create Greens functions
//============================================================================//

//============================================================================//
// Unbounded domain
//============================================================================//
	if(domain_bc[0] == 0 && domain_bc[1] == 0 && domain_bc[2] == 0)
	{
//----------------------------------------------------------------------------//
// Construct pencils
//----------------------------------------------------------------------------//
		class_pencil pen(true,false,false);

//----------------------------------------------------------------------------//
// FFT x-pencils
//----------------------------------------------------------------------------//
		ncell[0] = xpen_ext[rank].ncell[0];
		ncell[1] = xpen_ext[rank].ncell[1];
		ncell[2] = xpen_ext[rank].ncell[2];

		icell[0] = xpen_ext[rank].icell[0];
		icell[1] = xpen_ext[rank].icell[1];
		icell[2] = xpen_ext[rank].icell[2];

		dx[0] = xpen_ext[rank].dx[0];
		dx[1] = xpen_ext[rank].dx[1];
		dx[2] = xpen_ext[rank].dx[2];

		pen.resize( ncell[0] );

		mapG = new std::complex<double>[ncell[0]*ncell[1]*ncell[2]];

		xmin[0] = - ( 0.5 * double( xpen_ext[rank].ncell[0] ) )*dx[0];
		xmin[1] = - ( 0.5 * double( ypen_ext[rank].ncell[1] ) )*dx[1];
		xmin[2] = - ( 0.5 * double( ypen_ext[rank].ncell[2] ) )*dx[2];


		sigma = dx[0]/pi; // Spectral
//		sigma = 2.0*dx[0]; // super-Gaussian

		for (k = 0; k < ncell[2]; ++k )
		{
			kn = k * ncell[1];
			z  = xmin[2] + double( icell[2] + k )*dx[2];
			for (j = 0; j < ncell[1]; ++j )
			{
				kjn = (kn + j) * ncell[0];
				y   = xmin[1] + double( icell[1] + j )*dx[1];
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store mesh array in x-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				for(i = 0; i < ncell[0]; ++i )
				{
					x = xmin[0] + double( icell[0] + i )*dx[0];

					r = sqrt(x*x + y*y + z*z);

// 10th order super-Gaussian
/*
					rho = r/sigma;
					if(r > 0.25*dx[0])
					{
						pen.X[i] = {0.0 * dx[0]*dx[1]*dx[2], 0.0};
					}
					else
					{
						pen.X[i] = {0.0 * dx[0]*dx[1]*dx[2], 0.0};
					}
*/

// Spectral

					rho = r/sigma;
					if(r > 0.25*dx[0])
					{
						pen.X[i] = {c_1_2pi2 * sine_int(rho)/r * dx[0]*dx[1]*dx[2], 0.0};
					}
					else
					{
						pen.X[i] = {c_1_2pi2/sigma * dx[0]*dx[1]*dx[2], 0.0};
					}

/*
					if(r > 0.25*dx[0])
					{
						pen.X[i] = { 0.0, 0.0};
					}
					else
					{
						pen.X[i] = {1.0, 0.0};
					}
*/

				}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// FFT y-pencils
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				pen.fft_shift( );
				pen.fft( );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store array
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				for (i = 0; i < ncell[0]; ++i )
				{
					ijk = kjn + i;
					mapG[ijk] = pen.X[i];
				}

			}
		}

//----------------------------------------------------------------------------//
// Map to y-pencils
//----------------------------------------------------------------------------//
		map( comm );

//----------------------------------------------------------------------------//
// FFT y-pencils
//----------------------------------------------------------------------------//
		ncell[0] = ypen_ext[rank].ncell[0];
		ncell[1] = ypen_ext[rank].ncell[1];
		ncell[2] = ypen_ext[rank].ncell[2];

		rhsG = new std::complex<double>[ncell[0]*ncell[1]*ncell[2]];
		pen.resize( ncell[1] );

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
					pen.X[j] = mapG[ijk];
				}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// FFT y-pencils
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				pen.fft_shift( );
				pen.fft( );

//----------------------------------------------------------------------------//
// Store Fourier space Greens function
//----------------------------------------------------------------------------//
				for (j = 0; j < ncell[1]; ++j )
				{
					ijk = (kn + j) * ncell[0] + i;
					mapG[ijk] = pen.X[j];
				}

			}
		}


//----------------------------------------------------------------------------//
// FFT z-pencils
//----------------------------------------------------------------------------//
		pen.resize( ncell[2] );
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
					pen.X[k] = mapG[ijk];
				}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// FFT z-pencils
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
				pen.fft_shift( );
				pen.fft( );

//----------------------------------------------------------------------------//
// Store Fourier space Greens function
//----------------------------------------------------------------------------//
				for(k = 0; k < ncell[2]; ++k )
				{
					ijk = (k * ncell[1] + j) * ncell[0] + i;
					rhsG[ijk] = pen.X[k];
				}
			}
		}

//----------------------------------------------------------------------------//
// De-allocate mapG
//----------------------------------------------------------------------------//
		delete [] mapG;
		mapG = NULL;

	}
//============================================================================//
// Periodic-periodic-unbounded domain
//============================================================================//
/*
	else if(domain_bc[0] == 1 && domain_bc[1] == 1 && domain_bc[2] == 0)
	{

//----------------------------------------------------------------------------//
// Calculate Greens function direcly in Fourier space on y-pencil partition
//----------------------------------------------------------------------------//
		ncell[0] = ypen_ext[rank].ncell[0];
		ncell[1] = ypen_ext[rank].ncell[1];
		ncell[2] = ypen_ext[rank].ncell[2];

		icell[0] = ypen_ext[rank].icell[0];
		icell[1] = ypen_ext[rank].icell[1];
		icell[2] = ypen_ext[rank].icell[2];

		dx[0]    = ypen_ext[rank].dx[0];
		dx[1]    = ypen_ext[rank].dx[1];
		dx[2]    = ypen_ext[rank].dx[2];

//----------------------------------------------------------------------------//
// Calculate 1D free-space Greens function and Fourier transform it
//----------------------------------------------------------------------------//
		class_pencil pen(true,false,false);

		pen.resize( ncell[2] );

		sigma = dx[2]/pi; // Spectral
		C = sigma/pi;
		for(k = 0; k < ncell[2]; ++k )
		{
			z = xmin[2] + double( k )*dx[2];

// Spectral
			rho = std::abs(z)/sigma;
			pen.X[j] = {- C * ( sine_int( rho )*rho + cos(rho) ) * dx[2], 0.0};
		}

		pen.fft_shift( );
		pen.fft( );

//----------------------------------------------------------------------------//
// Construct Greens function in Fourier space on y-pencil partition
//----------------------------------------------------------------------------//
		rhsG = new std::complex<double>[ncell[0]*ncell[1]*ncell[2]];

		for (k = 0; k < ncell[2]; ++k )
		{
			kn = k * ncell[1];
			for (j = 0; j < ncell[1]; ++j )
			{
				kjn = (kn + j) * ncell[0];
				for(i = 0; i < ncell[0]; ++i )
				{
					ijk = kjn + i;

					if( icell[0] == 0 && i == 0 && j == 0 )
					{
						rhsG[ijk] = pen.X[k];
					}
					else
					{
						rhsG[ijk] = -1.0/( ikX[i]*ikX[i] + ikY[j]*ikY[j] + ikZ[k]*ikZ[k] );
					}

				}
			}
		}

	}
*/
//============================================================================//
// Periodic domain
//============================================================================//
	else
	{

//----------------------------------------------------------------------------//
// Calculate Greens function direcly in Fourier space on y-pencil partition
//----------------------------------------------------------------------------//
		ncell[0] = ypen_ext[rank].ncell[0];
		ncell[1] = ypen_ext[rank].ncell[1];
		ncell[2] = ypen_ext[rank].ncell[2];

		icell[0] = ypen_ext[rank].icell[0];
		icell[1] = ypen_ext[rank].icell[1];
		icell[2] = ypen_ext[rank].icell[2];

		dx[0]    = ypen_ext[rank].dx[0];
		dx[1]    = ypen_ext[rank].dx[1];
		dx[2]    = ypen_ext[rank].dx[2];

		rhsG = new std::complex<double>[ncell[0]*ncell[1]*ncell[2]];

		for (k = 0; k < ncell[2]; ++k )
		{
			kn = k * ncell[1];
			for (j = 0; j < ncell[1]; ++j )
			{
				kjn = (kn + j) * ncell[0];
				for(i = 0; i < ncell[0]; ++i )
				{
					ijk = kjn + i;

					if( icell[0] == 0 && i == 0 && j == 0  && k == 0 )
					{
						rhsG[ijk] = {0.0,0.0};
					}
					else
					{
						rhsG[ijk] = -1.0/( ikX[i]*ikX[i] + ikY[j]*ikY[j] + ikZ[k]*ikZ[k] );
					}

				}
			}
		}

	}

/*
	else
	{
		std::cerr << " [greenfish.solve]: Boundary condition configuration unknown."
		          << std::endl;
	}
*/

//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}


