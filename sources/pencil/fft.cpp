//----------------------------------------------------------------------------//
/*
  File:         fft.cpp

  Description:  Glassman's general N FFT algorithm see W. E. Ferguson: 
                "A simple derivation of Glassman's general N fast Fourier 
                transform", Comp. & Math. With Appls. 1982
*/
//----------------------------------------------------------------------------//

void class_pencil::fft( void )
{
//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//
	const double pi = acos(-1.0);

//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int i,j;
	int a, b, c;
	int ia, ib, ic, jc;
	bool inu;

	double angle;
	std::complex<double> delta, omega;
	std::complex<double> sumX,sumY,sumZ;

//----------------------------------------------------------------------------//
// Glassman FFT algorithm
//----------------------------------------------------------------------------//
	a = 1;
	b = nfft;
	c = 1;
	inu = true;

	CONTINUE_10:

	if(b > 1){ goto CONTINUE_30; }

	if(inu){ goto CONTINUE_END; }

	for (i = 0; i < nfft ; ++i )
	{
		if( bX )
		{
			X[i] = auxX[i];
		}
		if( bY )
		{
			Y[i] = auxY[i];
		}
		if( bZ )
		{
			Z[i] = auxZ[i];
		}
	}
	goto CONTINUE_END;

	CONTINUE_30:

	a = c*a;

	for (c = 2; c <= b ; ++c )
	{
		if( b % c == 0){ goto CONTINUE_50; }
	}

	CONTINUE_50:

	b = b/c;

	angle = 2.0*pi/double(a*c);

	delta = {cos(angle),-sin(angle)};
	omega = {1.0,0.0};


	if(inu)
	{ 
		for ( ic = 1; ic <= c ; ++ic )
		{
			for ( ia = 1; ia <= a ; ++ia )
			{
				for ( ib = 1; ib <= b ; ++ib )
				{

					if( bX )
					{
						sumX = X[ib + b*(c - 1) + b*c*(ia - 1) - 1];
					}
					if( bY )
					{
						sumY = Y[ib + b*(c - 1) + b*c*(ia - 1) - 1];
					}
					if( bZ )
					{
						sumZ = Z[ib + b*(c - 1) + b*c*(ia - 1) - 1];
					}

					for ( j = 2; j <= c ; ++j )
					{
						jc = c + 1 - j;
						if( bX )
						{
							sumX = X[ib + b*(jc-1) + b*c*(ia-1) - 1] + omega * sumX;
						}
						if( bY )
						{
							sumY = Y[ib + b*(jc-1) + b*c*(ia-1) - 1] + omega * sumY;
						}
						if( bZ )
						{
							sumZ = Z[ib + b*(jc-1) + b*c*(ia-1) - 1] + omega * sumZ;
						}
					}

					if( bX )
					{
						auxX[ib + b*(ia - 1) + b*a*(ic - 1) - 1] = sumX;
					}
					if( bY )
					{
						auxY[ib + b*(ia - 1) + b*a*(ic - 1) - 1] = sumY;
					}
					if( bZ )
					{
						auxZ[ib + b*(ia - 1) + b*a*(ic - 1) - 1] = sumZ;
					}

				}
				omega = delta * omega;
			}
		}
	}
	else
	{
		for ( ic = 1; ic <= c ; ++ic )
		{
			for ( ia = 1; ia <= a ; ++ia )
			{
				for ( ib = 1; ib <= b ; ++ib )
				{

					if( bX )
					{
						sumX = auxX[ib + b*(c - 1) + b*c*(ia - 1) - 1];
					}
					if( bY )
					{
						sumY = auxY[ib + b*(c - 1) + b*c*(ia - 1) - 1];
					}
					if( bZ )
					{
						sumZ = auxZ[ib + b*(c - 1) + b*c*(ia - 1) - 1];
					}

					for ( j = 2; j <= c ; ++j )
					{
						jc = c + 1 - j;
						if( bX )
						{
							sumX = auxX[ib + b*(jc-1) + b*c*(ia-1) - 1] + omega * sumX;
						}
						if( bY )
						{
							sumY = auxY[ib + b*(jc-1) + b*c*(ia-1) - 1] + omega * sumY;
						}
						if( bZ )
						{
							sumZ = auxZ[ib + b*(jc-1) + b*c*(ia-1) - 1] + omega * sumZ;
						}
					}

					if( bX )
					{
						X[ib + b*(ia - 1) + b*a*(ic - 1) - 1] = sumX;
					}
					if( bY )
					{
						Y[ib + b*(ia - 1) + b*a*(ic - 1) - 1] = sumY;
					}
					if( bZ )
					{
						Z[ib + b*(ia - 1) + b*a*(ic - 1) - 1] = sumZ;
					}

				}
				omega = delta * omega;
			}
		}
	}

	inu = !inu;

	goto CONTINUE_10;

	CONTINUE_END:

//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}
