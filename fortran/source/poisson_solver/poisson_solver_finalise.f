!------------------------------------------------------------------------------!
!  
!  File:         poisson_solver_finalise.f
!  
!  Description:  Finalises the poisson_solver
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_finalise( ps )


IMPLICIT NONE

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps
!	WRITE(*,*) " [poisson_solver]: Destructing Poisson solver object."

!------------------------------------------------------------------------------!
! De-allocate arrays from memory
!------------------------------------------------------------------------------!
	IF ( ALLOCATED(poisson_solver(ps)%rhsX) ) THEN
		DEALLOCATE(poisson_solver(ps)%rhsX)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%rhsY) ) THEN
		DEALLOCATE(poisson_solver(ps)%rhsY)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%rhsZ) ) THEN
		DEALLOCATE(poisson_solver(ps)%rhsZ)
	END IF

	IF ( ALLOCATED(poisson_solver(ps)%lhsX) ) THEN
		DEALLOCATE(poisson_solver(ps)%lhsX)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%lhsY) ) THEN
		DEALLOCATE(poisson_solver(ps)%lhsY)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%lhsZ) ) THEN
		DEALLOCATE(poisson_solver(ps)%lhsZ)
	END IF

	IF ( ALLOCATED(poisson_solver(ps)%ikX) ) THEN
		DEALLOCATE(poisson_solver(ps)%ikX)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%ikY) ) THEN
		DEALLOCATE(poisson_solver(ps)%ikY)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%ikZ) ) THEN
		DEALLOCATE(poisson_solver(ps)%ikZ)
	END IF

	IF ( ALLOCATED(poisson_solver(ps)%zeta) ) THEN
		DEALLOCATE(poisson_solver(ps)%zeta)
	END IF

	IF ( ALLOCATED(poisson_solver(ps)%mapG) ) THEN
		DEALLOCATE(poisson_solver(ps)%mapG)
	END IF

	IF ( ALLOCATED(poisson_solver(ps)%G2D) ) THEN
		DEALLOCATE(poisson_solver(ps)%G2D)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%G3D) ) THEN
		DEALLOCATE(poisson_solver(ps)%G3D)
	END IF

	IF ( ALLOCATED(poisson_solver(ps)%partition) ) THEN
		DEALLOCATE(poisson_solver(ps)%partition)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%xpen) ) THEN
		DEALLOCATE(poisson_solver(ps)%xpen)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%ypen) ) THEN
		DEALLOCATE(poisson_solver(ps)%ypen)
	END IF
	IF ( ALLOCATED(poisson_solver(ps)%zpen) ) THEN
		DEALLOCATE(poisson_solver(ps)%zpen)
	END IF

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE poisson_solver_finalise
