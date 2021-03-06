!------------------------------------------------------------------------------!
!  
!  File:         poisson_solver_initialise.f
!  
!  Description:  Finalises the poisson_solver
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_initialise( nsolvers )


IMPLICIT NONE

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: nsolvers

!------------------------------------------------------------------------------!
! Allocate solvers
!------------------------------------------------------------------------------!
	ALLOCATE( poisson_solver(nsolvers) )

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE poisson_solver_initialise
