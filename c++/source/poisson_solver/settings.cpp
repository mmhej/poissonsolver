//----------------------------------------------------------------------------//
/*
  File:         settings.cpp

  Description:  Provides set routines
*/
//----------------------------------------------------------------------------//

void poisson_solver::set_return_grad( bool toggle )
{
	lhs_grad = toggle;
	return;
}

void poisson_solver::set_return_div( bool toggle )
{
	lhs_div = toggle;
	return;
}

void poisson_solver::set_return_curl( bool toggle )
{
	lhs_curl = toggle;
	return;
}
void poisson_solver::set_reprojection( bool toggle )
{
	rhs_reproject = toggle;
	return;
}

