//----------------------------------------------------------------------------//
/*
  File:         destructor.cpp

  Description:  Destructor for class_pencil
*/
//----------------------------------------------------------------------------//
class_pencil::~class_pencil( void )
{

//----------------------------------------------------------------------------//
// De-allocate arrays from memory
//----------------------------------------------------------------------------//
	if( X != NULL )
	{
		delete [] X;
		X = NULL;
	}
	if( Y != NULL )
	{
		delete [] Y;
		Y = NULL;
	}
	if( Z != NULL )
	{
		delete [] Z;
		Z = NULL;
	}

	if( auxX != NULL )
	{
		delete [] auxX;
		auxX = NULL;
	}
	if( auxY != NULL )
	{
		delete [] auxY;
		auxY = NULL;
	}
	if( auxZ != NULL )
	{
		delete [] auxZ;
		auxZ = NULL;
	}

//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}
