!------------------------------------------------------------------------------!
!  
!  File:         resize.cpp
!  
!  Description:  Resizes the pencil arrays
!  
!------------------------------------------------------------------------------!
SUBROUTINE pencil_resize( pen, n )

IMPLICIT NONE

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	TYPE(class_pencil), INTENT(INOUT) :: pen
	INTEGER, INTENT(IN) :: n

!------------------------------------------------------------------------------!
! Set fft size
!------------------------------------------------------------------------------!
	pen%nfft = n;

!------------------------------------------------------------------------------!
! De-allocate fft vectors
!------------------------------------------------------------------------------!
	IF ( ALLOCATED(pen%X) ) THEN
		DEALLOCATE( pen%X )
	END IF
	IF ( ALLOCATED(pen%auxX) ) THEN
		DEALLOCATE( pen%auxX )
	END IF
	IF ( ALLOCATED(pen%Y) ) THEN
		DEALLOCATE( pen%Y )
	END IF
	IF ( ALLOCATED(pen%auxY) ) THEN
		DEALLOCATE( pen%auxY )
	END IF
	IF ( ALLOCATED(pen%Z) ) THEN
		DEALLOCATE( pen%Z )
	END IF
	IF ( ALLOCATED(pen%auxZ) ) THEN
		DEALLOCATE( pen%auxZ )
	END IF

!------------------------------------------------------------------------------!
! Re-allocate fft vectors
!------------------------------------------------------------------------------!
	IF ( pen%bX ) THEN
		ALLOCATE( pen%X(0:n-1), pen%auxX(0:n-1) )
		pen%X = CMPLX(0.0_MK,0.0_MK,MKC)
		pen%auxX = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( pen%bY ) THEN
		ALLOCATE( pen%Y(0:n-1), pen%auxY(0:n-1) )
		pen%Y = CMPLX(0.0_MK,0.0_MK,MKC)
		pen%auxY = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( pen%bZ ) THEN
		ALLOCATE( pen%Z(0:n-1), pen%auxZ(0:n-1) )
		pen%Z = CMPLX(0.0_MK,0.0_MK,MKC)
		pen%auxZ = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
9999 CONTINUE

RETURN
END SUBROUTINE pencil_resize
