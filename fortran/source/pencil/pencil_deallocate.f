!------------------------------------------------------------------------------!
!  
!  File:         pencil_deallocate.cpp
!  
!  Description:  Destructor for class_pencil
!  
!------------------------------------------------------------------------------!
SUBROUTINE pencil_deallocate( pen )

IMPLICIT NONE

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	TYPE(class_pencil), INTENT(INOUT) :: pen

!------------------------------------------------------------------------------!
! De-allocate arrays from memory
!------------------------------------------------------------------------------!
	IF ( ALLOCATED(pen%X) ) THEN
		DEALLOCATE( pen%X )
	END IF
	IF ( ALLOCATED(pen%Y) ) THEN
		DEALLOCATE( pen%Y )
	END IF
	IF ( ALLOCATED(pen%Z) ) THEN
		DEALLOCATE( pen%Z )
	END IF

	IF ( ALLOCATED(pen%auxX) ) THEN
		DEALLOCATE( pen%auxX )
	END IF
	IF ( ALLOCATED(pen%auxY) ) THEN
		DEALLOCATE( pen%auxY )
	END IF
	IF ( ALLOCATED(pen%auxZ) ) THEN
		DEALLOCATE( pen%auxZ )
	END IF

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
9999 CONTINUE

RETURN
END SUBROUTINE pencil_deallocate
