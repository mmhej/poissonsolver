//----------------------------------------------------------------------------//
/*
  File:         destructor.cpp

  Description:  
*/
//----------------------------------------------------------------------------//
class_greenfish::~class_greenfish( void )
{

//	std::cout << " [greenfish]: Destructing greenfish object." << std::endl;

//----------------------------------------------------------------------------//
// De-allocate arrays from memory
//----------------------------------------------------------------------------//
	if( rhsX != NULL )
	{
		delete [] rhsX;
		rhsX = NULL;
	}
	if( rhsY != NULL )
	{
		delete [] rhsY;
		rhsY = NULL;
	}
	if( rhsZ != NULL )
	{
		delete [] rhsZ;
		rhsZ = NULL;
	}

	if( lhsX != NULL )
	{
		delete [] lhsX;
		lhsX = NULL;
	}
	if( lhsY != NULL )
	{
		delete [] lhsY;
		lhsY = NULL;
	}
	if( lhsZ != NULL )
	{
		delete [] lhsZ;
		lhsZ = NULL;
	}

	if( ikX != NULL )
	{
		delete [] ikX;
		ikX = NULL;
	}
	if( ikY != NULL )
	{
		delete [] ikY;
		ikY = NULL;
	}
	if( ikZ != NULL )
	{
		delete [] ikZ;
		ikZ = NULL;
	}

	if( zeta != NULL )
	{
		delete [] zeta;
		zeta = NULL;
	}

	if( mapG != NULL )
	{
		delete [] mapG;
		mapG = NULL;
	}
	if( G2D != NULL )
	{
		delete [] G2D;
		G2D = NULL;
	}
	if( G3D != NULL )
	{
		delete [] G3D;
		G3D = NULL;
	}

//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}
