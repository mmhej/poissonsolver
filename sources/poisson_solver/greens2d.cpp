//----------------------------------------------------------------------------//
/*
  File:         greens2d.cpp

  Description:  Calculates the 2D Greens function and other convolution kernels
*/
//----------------------------------------------------------------------------//
void poisson_solver::greens2d( int dom_ncell[3], int dom_bc[3],
                               double dom_dx[3], int reg_order )
{

//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//
	const double pi = acos(-1.0);
	const double gamma = 0.5772156649015329;

	const double c_1_2pi = 0.5/pi;

//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int     nproc, rank;
	int     i, j, jn, ij, n;
	int     nfft;
	int     ncell[2];
	int     icell[2];
	double  dx[2];
	double  xmin[2];

	double  x,y,r,rho;
	double  arg, sum;
	double  sigma;
	double  h, dk;

	double * c;
	double  C1, C2;

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

	ncell[0] = dom_ncell[0];
	ncell[1] = dom_ncell[1];

	dx[0] = dom_dx[0];
	dx[1] = dom_dx[1];

	h = std::max(dx[0],dx[1]);

	xpen_ext = partition_setup( 0, dom_ncell, dom_bc, dom_dx, true, true, true );
	ypen_ext = partition_setup( 1, dom_ncell, dom_bc, dom_dx, true, true, true );

	class_communication comm( xpen_ext, ypen_ext );

//============================================================================//
// Create spectral derivative operator on y-pencil partition
//============================================================================//
	ncell[0] = ypen_ext[rank].ncell[0];
	ncell[1] = ypen_ext[rank].ncell[1];

	icell[0] = ypen_ext[rank].icell[0];
	icell[1] = ypen_ext[rank].icell[1];

	nfft = xpen_ext[rank].ncell[0];

	dk = 1.0/( double( nfft ) * dx[0] );
	ikX = new std::complex<double>[ncell[0]]();

	for (i = 0; i < ncell[0]; ++i )
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

	dk = 1.0/( double(ncell[1]) * dx[1] );
	ikY = new std::complex<double>[ncell[1]]();

	for (j = 0; j < ncell[1]; ++j )
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

//============================================================================//
// Create regularisation kernel
//============================================================================//
	sigma = 2.0*h;
	if( reg_order > 0 )
	{
		if( reg_order % 2 != 0 ){ reg_order += 1; }

		zeta = new double[ncell[0]*ncell[1]];

		c = new double[reg_order/2 - 1];
		c[0] = 1.0;
		for( n = 1; n <= (reg_order/2 - 1); ++n)
		{
			c[n] = c[n-1] * 1.0/double( n );
		}

		for(j = 0; j < ncell[1]; ++j )
		{
			jn = j * ncell[0];
			for (i = 0; i < ncell[0]; ++i )
			{
				ij = jn + i;
				arg = 0.5 * pow( sigma, 2 ) * ( pow( std::imag(ikX[i]), 2 )
				                              + pow( std::imag(ikY[j]), 2 ) );
				sum = 0.0;
				for( n = 0; n <= (reg_order/2 - 1); ++n)
				{
					sum += c[n] * pow(arg,n);
				}
				zeta[ij] = sum * exp( -arg );
			}
		}
	}

//============================================================================//
// Create Greens functions
//============================================================================//
//= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
// Unbounded domain
//= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	if(dom_bc[0] == 0 && dom_bc[1] == 0)
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

		icell[0] = xpen_ext[rank].icell[0];
		icell[1] = xpen_ext[rank].icell[1];

		dx[0] = xpen_ext[rank].dx[0];
		dx[1] = xpen_ext[rank].dx[1];

		pen.resize( ncell[0] );

		mapG = new std::complex<double>[ncell[0]*ncell[1]];

		xmin[0] = - ( 0.5 * double( xpen_ext[rank].ncell[0] ) )*dx[0];
		xmin[1] = - ( 0.5 * double( ypen_ext[rank].ncell[1] ) )*dx[1];

		sigma = h/pi; // Spectral
		C2    = log(2.0*sigma) - gamma;

		for (j = 0; j < ncell[1]; ++j )
		{
			jn = j * ncell[0];
			y  = xmin[1] + double( icell[1] + j )*dx[1];

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Store fft mesh array in x-pencil
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
			for(i = 0; i < ncell[0]; ++i )
			{
				x = xmin[0] + double( icell[0] + i )*dx[0];
				r = sqrt(x*x + y*y);

// Greens function
				rho = r/sigma;
				pen.X[i] = {- c_1_2pi * ( bessel_int_J0( rho ) + C2 ) * dx[0]*dx[1], 0.0};
			}

			pen.fft_shift( );
			pen.fft( );

			for (i = 0; i < ncell[0]; ++i )
			{
				ij = jn + i;
				mapG[ij] = pen.X[i];
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

		G2D = new std::complex<double>[ncell[0]*ncell[1]];
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
				G2D[ij] = pen.X[j];
			}

		}

//----------------------------------------------------------------------------//
// De-allocate mapG
//----------------------------------------------------------------------------//
		delete [] mapG;
		mapG = NULL;

	}
//= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
// Periodic domain
//= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	else if(dom_bc[0] == 1 && dom_bc[1] == 1)
	{

//----------------------------------------------------------------------------//
// Calculate Greens function direcly in Fourier space on y-pencil partition
//----------------------------------------------------------------------------//
		ncell[0] = ypen_ext[rank].ncell[0];
		ncell[1] = ypen_ext[rank].ncell[1];

		icell[0] = ypen_ext[rank].icell[0];
		icell[1] = ypen_ext[rank].icell[1];

		dx[0]    = ypen_ext[rank].dx[0];
		dx[1]    = ypen_ext[rank].dx[1];

		G2D = new std::complex<double>[ncell[0]*ncell[1]];

		for(j = 0; j < ncell[1]; ++j )
		{
			jn = j * ncell[0];
			for (i = 0; i < ncell[0]; ++i )
			{
				ij = jn + i;
				if( icell[0] == 0 && i == 0 && j == 0 )
				{
					G2D[ij] = {0.0,0.0};
				}
				else
				{
					G2D[ij] = -1.0/( ikX[i]*ikX[i] + ikY[j]*ikY[j] );
				}
			}
		}


	}
//= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
// Periodic-unbounded domain
//= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	else if(dom_bc[0] == 1 && dom_bc[1] == 0)
	{

		ncell[0] = ypen_ext[rank].ncell[0];
		ncell[1] = ypen_ext[rank].ncell[1];

		icell[0] = ypen_ext[rank].icell[0];
		icell[1] = ypen_ext[rank].icell[1];

		dx[0]    = ypen_ext[rank].dx[0];
		dx[1]    = ypen_ext[rank].dx[1];

		xmin[0] = - ( 0.5 * double( xpen_ext[rank].ncell[0] ) )*dx[0];
		xmin[1] = - ( 0.5 * double( ypen_ext[rank].ncell[1] ) )*dx[1];

//----------------------------------------------------------------------------//
// Calculate 1D free-space Greens function and Fourier transform it
//----------------------------------------------------------------------------//
		class_pencil pen(true,false,false);

		pen.resize( ncell[1] );

		sigma = dx[1]/pi; // Spectral
		C1 = sigma/pi;
		for(j = 0; j < ncell[1]; ++j )
		{
			y = xmin[1] + double( j )*dx[1];

// Spectral
			rho = std::abs(y)/sigma;
			pen.X[j] = {- C1 * ( sine_int( rho )*rho + cos(rho) ) * dx[1], 0.0};
		}

		pen.fft_shift( );
		pen.fft( );

//----------------------------------------------------------------------------//
// Construct Greens function in Fourier space on y-pencil partition
//----------------------------------------------------------------------------//
		G2D = new std::complex<double>[ncell[0]*ncell[1]];
		for(j = 0; j < ncell[1]; ++j )
		{
			jn = j * ncell[0];
			for (i = 0; i < ncell[0]; ++i )
			{
				ij = jn + i;

				if( icell[0] == 0 && i == 0 )
				{
					G2D[ij] = pen.X[j];
				}
				else
				{
					G2D[ij] = -1.0/( ikX[i]*ikX[i] + ikY[j]*ikY[j] );
				}

			}
		}

	}
//= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
// Unbounded-periodic domain
//= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	else if(dom_bc[0] == 0 && dom_bc[1] == 1)
	{

		ncell[0] = ypen_ext[rank].ncell[0];
		ncell[1] = ypen_ext[rank].ncell[1];

		icell[0] = ypen_ext[rank].icell[0];
		icell[1] = ypen_ext[rank].icell[1];

		dx[0]    = ypen_ext[rank].dx[0];
		dx[1]    = ypen_ext[rank].dx[1];

		nfft = xpen_ext[rank].ncell[0];
		xmin[0] = - ( 0.5 * double( nfft ) )*dx[0];

//----------------------------------------------------------------------------//
// Calculate 1D free-space Greens function and Fourier transform it
//----------------------------------------------------------------------------//
		class_pencil pen(true,false,false);

		pen.resize( nfft );

		sigma = dx[0]/pi; // Spectral
		C1 = sigma/pi;
		for(i = 0; i < nfft; ++i )
		{
			x = xmin[0] + double( i )*dx[0];

// Spectral
			rho = std::abs(x)/sigma;
			pen.X[j] = {- C1 * ( sine_int( rho )*rho + cos(rho) ) * dx[0], 0.0};
		}

		pen.fft_shift( );
		pen.fft( );

//----------------------------------------------------------------------------//
// Construct Greens function in Fourier space on y-pencil partition
//----------------------------------------------------------------------------//
		G2D = new std::complex<double>[ncell[0]*ncell[1]];
		for(j = 0; j < ncell[1]; ++j )
		{
			jn = j * ncell[0];
			for (i = 0; i < ncell[0]; ++i )
			{
				ij = jn + i;

				if( j == 0 )
				{
					G2D[ij] = pen.X[icell[0] + i];
				}
				else
				{
					G2D[ij] = -1.0/( ikX[i]*ikX[i] + ikY[j]*ikY[j] );
				}

			}
		}

	}
	else
	{
		std::cerr << " [poisson_solver.greens2d]: Boundary condition configuration unknown."
		          << std::endl;
	}


//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}


