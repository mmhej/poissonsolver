//----------------------------------------------------------------------------//
/*
  File:         poisson_solver.cpp

  Description:  Includes the source code for poisson_solver
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Include header file
//----------------------------------------------------------------------------//
#include "poisson_solver.hpp"

//----------------------------------------------------------------------------//
// Include subroutines
//----------------------------------------------------------------------------//
#include "poisson_solver/constructor.cpp"
#include "poisson_solver/destructor.cpp"

#include "poisson_solver/setup2d.cpp"
#include "poisson_solver/setup3d.cpp"

#include "poisson_solver/settings.cpp"

#include "poisson_solver/map.cpp"

#include "poisson_solver/special_functions.cpp"

#include "poisson_solver/greens2d.cpp"
#include "poisson_solver/greens3d.cpp"

#include "poisson_solver/solve2d.cpp"
#include "poisson_solver/solve3d.cpp"

#include "poisson_solver/push.cpp"
#include "poisson_solver/pull.cpp"

//----------------------------------------------------------------------------//
// Include customised subroutines
//----------------------------------------------------------------------------//
//#include "../custom/custom_solve2d.cpp"


