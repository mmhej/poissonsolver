//----------------------------------------------------------------------------//
/*
  File:         resize.cpp

  Description:  Resizes the pencil arrays
*/
//----------------------------------------------------------------------------//

void class_pencil::resize( int n )
{
//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Objects
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Set fft size
//----------------------------------------------------------------------------//
	nfft = n;

//----------------------------------------------------------------------------//
// Re-size fft vectors
//----------------------------------------------------------------------------//
	if( bX )
	{
		X    = new std::complex<double>[n]();
		auxX = new std::complex<double>[n]();
	}
	if( bY )
	{
		Y    = new std::complex<double>[n]();
		auxY = new std::complex<double>[n]();
	}
	if( bZ )
	{
		Z    = new std::complex<double>[n]();
		auxZ = new std::complex<double>[n]();
	}



//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}