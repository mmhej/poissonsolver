//----------------------------------------------------------------------------//
/*
  File:         greens2d.cpp

  Description:  
*/
//----------------------------------------------------------------------------//
void class_greenfish::greens2d(  )
{

//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//
	const double pi = acos(-1.0);

	const double c_1_2pi = 0.5/pi;
	const double c_sqrt2 = sqrt(2.0);
	const double c1 =   25.0/24.0;
	const double c2 = - 23.0/48.0;
	const double c3 =   13.0/192.0;
	const double c4 = -  1.0/384.0;
	const double gamma = 0.5772156649015329;

//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int     nproc, rank;
	int     i, j, ij;
	int     nfft;
	int     ncell[2];
	int     icell[2];
	int     bound_cond[2];
	double  dx[2];
	double  xmin[2];

	double  x,y,r,rho;
	double sigma;
	double dk;
	double * kX;
	double * kY;

	double kL;
	double kx2, ky2;

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//============================================================================//
// Create extended partitions and communications
//============================================================================//
	class_communication comm;
	std::vector<class_partition_info> xpen;
	std::vector<class_partition_info> ypen;

	ncell[0] = domain_ncell[0];
	ncell[1] = domain_ncell[1];

	bound_cond[0] = domain_bc[0];
	bound_cond[1] = domain_bc[1];

	dx[0] = domain_dx[0];
	dx[1] = domain_dx[1];

	xpen = partition_setup2d( 0, ncell, bound_cond, dx, true, true );
	ypen = partition_setup2d( 1, ncell, bound_cond, dx, true, true );

	comm.setup2d( xpen, ypen );


//============================================================================//
// Create spectral derivative operator on y-pencil partition
//============================================================================//
	ncell[0] = ypen[rank].ncell[0];
	ncell[1] = ypen[rank].ncell[1];

	icell[0] = ypen[rank].icell[0];
	icell[1] = ypen[rank].icell[1];

	dk = 1.0/double( xpen[rank].ncell[0] );
	ikX = new std::complex<double>[ncell[0]]();

	for (i = 0; i < ncell[0]; ++i )
	{
		if( 2*i < ncell[0] )
		{
			ikX[i] = {0.0, 2.0 * pi * double(icell[0]+i) * dk};
		}
		else
		{
			ikX[i] = {0.0, 2.0 * pi * double(icell[0]+i) * dk - 1.0};
		}
	}

	dk = 1.0/double(ncell[1]);
	ikY = new std::complex<double>[ncell[1]]();

	for (j = 0; j < ncell[1]; ++j )
	{
		if( 2*j < ncell[1] )
		{
			ikY[j] = {0.0, 2.0 * pi * double(j) * dk};
		}
		else
		{
			ikY[j] = {0.0, 2.0 * pi * double(j) * dk - 1.0};
		}
	}


//============================================================================//
// Create Greens functions
//============================================================================//
	if(domain_bc[0] == 0 && domain_bc[1] == 0)
	{
//----------------------------------------------------------------------------//
// Construct pencils
//----------------------------------------------------------------------------//
		class_pencil pen(true,false,false);

//----------------------------------------------------------------------------//
// FFT x-pencils
//----------------------------------------------------------------------------//
		ncell[0] = xpen[rank].ncell[0];
		ncell[1] = xpen[rank].ncell[1];

		icell[0] = xpen[rank].icell[0];
		icell[1] = xpen[rank].icell[1];

		dx[0] = xpen[rank].dx[0];
		dx[1] = xpen[rank].dx[1];

		pen.resize( ncell[0] );

		mapG = new std::complex<double>[ncell[0]*ncell[1]];

		xmin[0] = - ( 0.5 * double( xpen[rank].ncell[0] ) )*dx[0];
		xmin[1] = - ( 0.5 * double( ypen[rank].ncell[1] ) )*dx[1];

		for (j = 0; j < ncell[1]; ++j )
		{
			y = xmin[1] + double( icell[1] + j )*dx[1];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft mesh array in x-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for(i = 0; i < ncell[0]; ++i )
			{
				x = xmin[0] + double( icell[0] + i )*dx[0];
				r = sqrt(x*x + y*y);


// 10th order super-gaussian
/*
				sigma = 2.0*dx[0];
				rho = r/sigma;
				if(r > 0.25*dx[0])
				{
					pen.X[i] = {- (log(r) + 0.5 * exp_int(1,0.5*pow(rho,2)) - (c1 + c2 * pow(rho,2) + c3 * pow(rho,4) + c4 * pow(rho,6) ) * exp(-0.5 * pow(rho,2)) ) * c_1_2pi * dx[0]*dx[1], 0.0};
				}
				else
				{
					pen.X[i] = {( 0.5*gamma - log(c_sqrt2 * sigma) + c1 )*c_1_2pi * dx[0]*dx[1], 0.0};
				}
*/

// Spectral
				sigma = dx[0]/pi;
				rho = r/sigma;
				pen.X[i] = {- c_1_2pi * ( bessel_int_J0( rho ) + log(2.0*sigma) - gamma ) * dx[0]*dx[1], 0.0};
			}

			pen.fft_shift( );
			pen.fft( );

			for (i = 0; i < ncell[0]; ++i )
			{
				ij = j * ncell[0] + i;
				mapG[ij] = pen.X[i];
			}

		}

//----------------------------------------------------------------------------//
// Map to y-pencils
//----------------------------------------------------------------------------//
		map2d( comm );

//----------------------------------------------------------------------------//
// FFT y-pencils
//----------------------------------------------------------------------------//
		ncell[0] = ypen[rank].ncell[0];
		ncell[1] = ypen[rank].ncell[1];

		rhsG = new std::complex<double>[ncell[0]*ncell[1]];
		pen.resize( ncell[1] );

		for (i = 0; i < ncell[0]; ++i )
		{

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft mesh array in y-pencil (fft-shifted)
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for(j = 0; j < ncell[1]; ++j )
			{
				ij = j * ncell[0] + i;
				pen.X[j] = mapG[ij];
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
				ij = j * ncell[0] + i;
				rhsG[ij] = pen.X[j];
			}

		}

//----------------------------------------------------------------------------//
// De-allocate mapG
//----------------------------------------------------------------------------//
		delete [] mapG;
		mapG = NULL;

	}
	else if(domain_bc[0] == 1 && domain_bc[1] == 1)
	{

//----------------------------------------------------------------------------//
// FFT y-pencils
//----------------------------------------------------------------------------//
		ncell[0] = ypen[rank].ncell[0];
		ncell[1] = ypen[rank].ncell[1];

		icell[0] = ypen[rank].icell[0];
		icell[1] = ypen[rank].icell[1];

		dx[0]    = ypen[rank].dx[0];
		dx[1]    = ypen[rank].dx[1];

		rhsG = new std::complex<double>[ncell[0]*ncell[1]];

		for (i = 0; i < ncell[0]; ++i )
		{
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft mesh array in y-pencil (fft-shifted)
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for(j = 0; j < ncell[1]; ++j )
			{
				ij = j * ncell[0] + i;

				if( j > 0 || i > 0 || icell[0] > 0 )
				{
					rhsG[ij] = -dx[0]*dx[1]/( ikX[i]*ikX[i] + ikY[j]*ikY[j] );
				}
				else
				{
					rhsG[ij] = {0.0,0.0};
				}
			}
		}

	}

//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}


