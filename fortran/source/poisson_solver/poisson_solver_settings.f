!------------------------------------------------------------------------------!
!  
!  File:         settings.f
!  Description:  Provides set routines
!  
!------------------------------------------------------------------------------!

SUBROUTINE poisson_solver_set_return_grad( toggle )

IMPLICIT NONE

	LOGICAL, INTENT(IN) :: toggle

	poisson_solver%lhs_grad = toggle

RETURN
END SUBROUTINE poisson_solver_set_return_grad



SUBROUTINE poisson_solver_set_return_div( toggle )

IMPLICIT NONE

	LOGICAL, INTENT(IN) :: toggle

	poisson_solver%lhs_div = toggle

RETURN
END SUBROUTINE poisson_solver_set_return_div



SUBROUTINE poisson_solver_set_return_curl( toggle )

IMPLICIT NONE

	LOGICAL, INTENT(IN) :: toggle

	poisson_solver%lhs_curl = toggle

RETURN
END SUBROUTINE poisson_solver_set_return_curl


SUBROUTINE poisson_solver_set_reprojection( toggle )

IMPLICIT NONE

	LOGICAL, INTENT(IN) :: toggle

	poisson_solver%rhs_reproject = toggle

RETURN
END SUBROUTINE poisson_solver_set_reprojection



