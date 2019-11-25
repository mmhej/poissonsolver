!------------------------------------------------------------------------------!
!  
!  File:         finalise.f
!  
!  Description:  Finalises the poisson_solver
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_finalise( )


IMPLICIT NONE

!	WRITE(*,*) " [poisson_solver]: Destructing Poisson solver object."

!------------------------------------------------------------------------------!
! De-allocate arrays from memory
!------------------------------------------------------------------------------!
	IF ( ALLOCATED(poisson_solver%rhsX) ) THEN
		DEALLOCATE(poisson_solver%rhsX)
	END IF
	IF ( ALLOCATED(poisson_solver%rhsY) ) THEN
		DEALLOCATE(poisson_solver%rhsY)
	END IF
	IF ( ALLOCATED(poisson_solver%rhsZ) ) THEN
		DEALLOCATE(poisson_solver%rhsZ)
	END IF

	IF ( ALLOCATED(poisson_solver%lhsX) ) THEN
		DEALLOCATE(poisson_solver%lhsX)
	END IF
	IF ( ALLOCATED(poisson_solver%lhsY) ) THEN
		DEALLOCATE(poisson_solver%lhsY)
	END IF
	IF ( ALLOCATED(poisson_solver%lhsZ) ) THEN
		DEALLOCATE(poisson_solver%lhsZ)
	END IF

	IF ( ALLOCATED(poisson_solver%ikX) ) THEN
		DEALLOCATE(poisson_solver%ikX)
	END IF
	IF ( ALLOCATED(poisson_solver%ikY) ) THEN
		DEALLOCATE(poisson_solver%ikY)
	END IF
	IF ( ALLOCATED(poisson_solver%ikZ) ) THEN
		DEALLOCATE(poisson_solver%ikZ)
	END IF

	IF ( ALLOCATED(poisson_solver%zeta) ) THEN
		DEALLOCATE(poisson_solver%zeta)
	END IF

	IF ( ALLOCATED(poisson_solver%mapG) ) THEN
		DEALLOCATE(poisson_solver%mapG)
	END IF

	IF ( ALLOCATED(poisson_solver%G2D) ) THEN
		DEALLOCATE(poisson_solver%G2D)
	END IF
	IF ( ALLOCATED(poisson_solver%G3D) ) THEN
		DEALLOCATE(poisson_solver%G3D)
	END IF

	IF ( ALLOCATED(poisson_solver%partition) ) THEN
		DEALLOCATE(poisson_solver%partition)
	END IF
	IF ( ALLOCATED(poisson_solver%xpen) ) THEN
		DEALLOCATE(poisson_solver%xpen)
	END IF
	IF ( ALLOCATED(poisson_solver%ypen) ) THEN
		DEALLOCATE(poisson_solver%ypen)
	END IF
	IF ( ALLOCATED(poisson_solver%zpen) ) THEN
		DEALLOCATE(poisson_solver%zpen)
	END IF


!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE poisson_solver_finalise
