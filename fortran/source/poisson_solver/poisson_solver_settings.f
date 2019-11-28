!------------------------------------------------------------------------------!
!  
!  File:         settings.f
!  Description:  Provides set routines
!  
!------------------------------------------------------------------------------!

SUBROUTINE poisson_solver_set_return_grad( ps, toggle )

IMPLICIT NONE

	INTEGER,  INTENT(IN) :: ps
	LOGICAL, INTENT(IN) :: toggle

	poisson_solver(ps)%lhs_grad = toggle

RETURN
END SUBROUTINE poisson_solver_set_return_grad



SUBROUTINE poisson_solver_set_return_div( ps, toggle )

IMPLICIT NONE

	INTEGER,  INTENT(IN) :: ps
	LOGICAL, INTENT(IN) :: toggle

	poisson_solver(ps)%lhs_div = toggle

RETURN
END SUBROUTINE poisson_solver_set_return_div



SUBROUTINE poisson_solver_set_return_curl( ps, toggle )

IMPLICIT NONE

	INTEGER,  INTENT(IN) :: ps
	LOGICAL, INTENT(IN) :: toggle

	poisson_solver(ps)%lhs_curl = toggle

RETURN
END SUBROUTINE poisson_solver_set_return_curl


SUBROUTINE poisson_solver_set_reprojection( ps, toggle )

IMPLICIT NONE

	INTEGER,  INTENT(IN) :: ps
	LOGICAL, INTENT(IN) :: toggle

	poisson_solver(ps)%rhs_reproject = toggle

RETURN
END SUBROUTINE poisson_solver_set_reprojection



