!------------------------------------------------------------------------------!
!  
!  File:         poisson_solver_solve2d.f
!  
!  Description:  Solves the Poisson equation in 2D
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_solve2d(ps)

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps
	REAL(MK), PARAMETER :: pi = 3.1415926535897932_MK

!------------------------------------------------------------------------------!
! Variables
!------------------------------------------------------------------------------!
	INTEGER :: nproc, rank, ierr
	INTEGER :: i, j, jn, ij
	INTEGER :: nfft
	INTEGER,DIMENSION(2) :: ncell

	LOGICAL :: rX,rY,rZ,lX,lY,lZ,mG

	TYPE(class_pencil) :: pen_rhs
	TYPE(class_pencil) :: pen_lhs

!------------------------------------------------------------------------------!
! Get MPI info
!------------------------------------------------------------------------------!
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

!------------------------------------------------------------------------------!
! Find out which fields to map
!------------------------------------------------------------------------------!
	rX = .FALSE.
	rY = .FALSE.
	rZ = .FALSE.
	lX = .FALSE.
	lY = .FALSE.
	lZ = .FALSE.

	IF ( ALLOCATED( poisson_solver(ps)%rhsX ) ) THEN
		rX = .TRUE.
	END IF
	IF ( ALLOCATED( poisson_solver(ps)%rhsY ) ) THEN
		rY = .TRUE.
	END IF
	IF ( ALLOCATED( poisson_solver(ps)%rhsZ ) ) THEN
		rZ = .TRUE.
	END IF

	IF (poisson_solver(ps)%lhs_grad) THEN
		IF (rX) THEN
			lX = .TRUE.
			lY = .TRUE.
		ELSE
			WRITE(*,'(2A)') " [poisson_solver_solve2d]: ", &
			     "Error rhs fields have not been pushed correctly for grad option." 
			GO TO 9999
		END IF
	ELSE IF (poisson_solver(ps)%lhs_div) THEN
		IF (rX .AND. rY) THEN
			lX = .TRUE.
		ELSE
			WRITE(*,'(2A)') " [poisson_solver_solve2d]: ", &
			     "Error rhs fields have not been pushed correctly for div option." 
			GO TO 9999
		END IF
	ELSE IF (poisson_solver(ps)%lhs_curl) THEN
		IF (rZ) THEN
			lX = .TRUE.
			lY = .TRUE.
		ELSE
			WRITE(*,'(2A)') " [poisson_solver_solve2d]: ", &
			     "Error rhs fields have not been pushed correctly for curl option." 
			GO TO 9999
		END IF
	ELSE
		IF (rX) THEN
			lX = .TRUE.
		END IF
		IF (rY) THEN
			lY = .TRUE.
		END IF
		IF (rZ) THEN
			lZ = .TRUE.
		END IF
	END IF

!------------------------------------------------------------------------------!
! Construct pencils
!------------------------------------------------------------------------------!
	CALL pencil_setup( pen_rhs, rX, rY, rZ)
	CALL pencil_setup( pen_lhs, lX, lY, lZ)

!------------------------------------------------------------------------------!
! FFT x-pencils
!------------------------------------------------------------------------------!
	ncell(1) = poisson_solver(ps)%xpen(rank)%ncell(1)
	ncell(2) = poisson_solver(ps)%xpen(rank)%ncell(2)

	CALL pencil_resize( pen_rhs, ncell(1) )

	DO j = 0,ncell(2)-1
		jn = j * ncell(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store fft mesh array in x-pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		DO i = 0,ncell(1)-1
			ij = jn + i
			IF ( rX ) THEN
				pen_rhs%X(i) = poisson_solver(ps)%rhsX(ij)
			END IF
			IF ( rY ) THEN
				pen_rhs%Y(i) = poisson_solver(ps)%rhsY(ij)
			END IF
			IF ( rZ ) THEN
				pen_rhs%Z(i) = poisson_solver(ps)%rhsZ(ij)
			END IF
		END DO

		CALL pencil_fft( pen_rhs )

		DO i = 0,ncell(1)-1
			ij = jn + i
			IF ( rX ) THEN
				poisson_solver(ps)%rhsX(ij) = pen_rhs%X(i)
			END IF
			IF ( rY ) THEN
				poisson_solver(ps)%rhsY(ij) = pen_rhs%Y(i)
			END IF
			IF ( rZ ) THEN
				poisson_solver(ps)%rhsZ(ij) = pen_rhs%Z(i)
			END IF
		END DO
	END DO

!------------------------------------------------------------------------------!
! Map to y-pencils
!------------------------------------------------------------------------------!
	CALL poisson_solver_map( ps, poisson_solver(ps)%xpen2ypen )

!------------------------------------------------------------------------------!
! FFT y-pencils and perform Fourier space operations
!------------------------------------------------------------------------------!
	ncell(1) = poisson_solver(ps)%ypen(rank)%ncell(1)
	ncell(2) = poisson_solver(ps)%ypen(rank)%ncell(2)

	IF (poisson_solver(ps)%bc(2) .EQ. 0) THEN
		nfft = 2*ncell(2)
	ELSE
		nfft = ncell(2)
	END IF

	CALL pencil_resize( pen_rhs, nfft )
	CALL pencil_resize( pen_lhs, nfft )

	IF ( lX ) THEN
		IF( ALLOCATED(poisson_solver(ps)%lhsX) ) THEN
			IF ( SIZE(poisson_solver(ps)%lhsX) .NE. ncell(1)*ncell(2) )THEN
				DEALLOCATE( poisson_solver(ps)%lhsX )
				ALLOCATE( poisson_solver(ps)%lhsX( 0:ncell(1)*ncell(2)-1 ) )
			END IF
		ELSE
			ALLOCATE( poisson_solver(ps)%lhsX( 0:ncell(1)*ncell(2)-1 ) )
		END IF
		poisson_solver(ps)%lhsX = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( lY ) THEN
		IF( ALLOCATED(poisson_solver(ps)%lhsY) ) THEN
			IF ( SIZE(poisson_solver(ps)%lhsY) .NE. ncell(1)*ncell(2) )THEN
				DEALLOCATE( poisson_solver(ps)%lhsY )
				ALLOCATE( poisson_solver(ps)%lhsY( 0:ncell(1)*ncell(2)-1 ) )
			END IF
		ELSE
			ALLOCATE( poisson_solver(ps)%lhsY( 0:ncell(1)*ncell(2)-1 ) )
		END IF
		poisson_solver(ps)%lhsY = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( lZ ) THEN
		IF( ALLOCATED(poisson_solver(ps)%lhsZ) ) THEN
			IF ( SIZE(poisson_solver(ps)%lhsZ) .NE. ncell(1)*ncell(2) )THEN
				DEALLOCATE( poisson_solver(ps)%lhsZ )
				ALLOCATE( poisson_solver(ps)%lhsZ( 0:ncell(1)*ncell(2)-1 ) )
			END IF
		ELSE
			ALLOCATE( poisson_solver(ps)%lhsZ( 0:ncell(1)*ncell(2)-1 ) )
		END IF
		poisson_solver(ps)%lhsZ = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF


	DO i = 0,ncell(1)-1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store fft mesh array in y-pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		DO j = 0,ncell(2)-1
			ij = j * ncell(1) + i
			IF ( rX ) THEN
				pen_rhs%X(j) = poisson_solver(ps)%rhsX(ij)
			END IF
			IF ( rY ) THEN
				pen_rhs%Y(j) = poisson_solver(ps)%rhsY(ij)
			END IF
			IF ( rZ ) THEN
				pen_rhs%Z(j) = poisson_solver(ps)%rhsZ(ij)
			END IF
		END DO
! Zero-padding
		IF ( nfft .GT. ncell(2) ) THEN
			DO j = ncell(2),nfft-1
				IF ( rX ) THEN
					pen_rhs%X(j) = CMPLX(0.0_MK,0.0_MK,MKC)
				END IF
				IF ( rY ) THEN
					pen_rhs%Y(j) = CMPLX(0.0_MK,0.0_MK,MKC)
				END IF
				IF ( rZ ) THEN
					pen_rhs%Z(j) = CMPLX(0.0_MK,0.0_MK,MKC)
				END IF
			END DO
		END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! FFT y-pencils
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		CALL pencil_fft( pen_rhs )

!------------------------------------------------------------------------------!
! Perform Fourier space operations
!------------------------------------------------------------------------------!
		DO j = 0,nfft-1
			ij = j * ncell(1) + i

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Regularise rhs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			IF ( poisson_solver(ps)%regularisation .GT. 0 ) THEN
				IF ( rX ) THEN
					pen_rhs%X(j) = pen_rhs%X(j)*poisson_solver(ps)%zeta(ij)
				END IF
				IF ( rY ) THEN
					pen_rhs%Y(j) = pen_rhs%Y(j)*poisson_solver(ps)%zeta(ij)
				END IF
				IF ( rZ ) THEN
					pen_rhs%Z(j) = pen_rhs%Z(j)*poisson_solver(ps)%zeta(ij)
				END IF
			END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Do convolution
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			IF (poisson_solver(ps)%lhs_grad) THEN
				pen_lhs%X(j) = poisson_solver(ps)%ikX(i) * poisson_solver(ps)%G2D(ij) * pen_rhs%X(j)
				pen_lhs%Y(j) = poisson_solver(ps)%ikY(j) * poisson_solver(ps)%G2D(ij) * pen_rhs%X(j)
			ELSE IF (poisson_solver(ps)%lhs_div) THEN
				pen_lhs%X(j) = poisson_solver(ps)%G2D(ij) * &
				                         ( poisson_solver(ps)%ikX(i) * pen_rhs%X(j) &
				                         + poisson_solver(ps)%ikY(j) * pen_rhs%Y(j) )
			ELSE IF (poisson_solver(ps)%lhs_curl) THEN
				pen_lhs%X(j) =   poisson_solver(ps)%ikY(j) * poisson_solver(ps)%G2D(ij) * pen_rhs%Z(j)
				pen_lhs%Y(j) = - poisson_solver(ps)%ikX(i) * poisson_solver(ps)%G2D(ij) * pen_rhs%Z(j)
			ELSE
				IF ( lX ) THEN
					pen_lhs%X(j) = poisson_solver(ps)%G2D(ij) * pen_rhs%X(j)
				END IF
				IF ( lY ) THEN
					pen_lhs%Y(j) = poisson_solver(ps)%G2D(ij) * pen_rhs%Y(j)
				END IF
				IF ( lZ ) THEN
					pen_lhs%Z(j) = poisson_solver(ps)%G2D(ij) * pen_rhs%Z(j)
				END IF
			END IF
		END DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! IFFT y-pencils
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		CALL pencil_ifft( pen_rhs )
		CALL pencil_ifft( pen_lhs )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store solution
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		DO j = 0,ncell(2)-1
			ij = j * ncell(1) + i
			IF ( rX ) THEN
				poisson_solver(ps)%rhsX(ij) = pen_rhs%X(j)/REAL(nfft,MK)
			END IF
			IF ( rY ) THEN
				poisson_solver(ps)%rhsY(ij) = pen_rhs%Y(j)/REAL(nfft,MK)
			END IF
			IF ( rZ ) THEN
				poisson_solver(ps)%rhsZ(ij) = pen_rhs%Z(j)/REAL(nfft,MK)
			END IF

			IF ( lX ) THEN
				poisson_solver(ps)%lhsX(ij) = pen_lhs%X(j)/REAL(nfft,MK)
			END IF
			IF ( lY ) THEN
				poisson_solver(ps)%lhsY(ij) = pen_lhs%Y(j)/REAL(nfft,MK)
			END IF
			IF ( lZ ) THEN
				poisson_solver(ps)%lhsZ(ij) = pen_lhs%Z(j)/REAL(nfft,MK)
			END IF
		END DO
	END DO

!------------------------------------------------------------------------------!
! Map to x-pencils
!------------------------------------------------------------------------------!
	CALL poisson_solver_map( ps, poisson_solver(ps)%ypen2xpen )

!------------------------------------------------------------------------------!
! IFFT x-pencils
!------------------------------------------------------------------------------!
	ncell(1) = poisson_solver(ps)%xpen(rank)%ncell(1)
	ncell(2) = poisson_solver(ps)%xpen(rank)%ncell(2)

	CALL pencil_resize( pen_rhs, ncell(1) )
	CALL pencil_resize( pen_lhs, ncell(1) )

	DO j = 0,ncell(2)-1
		jn = j * ncell(1)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store fft mesh array in x-pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		DO i = 0,ncell(1)-1
			ij = jn + i
			IF ( rX ) THEN
				pen_rhs%X(i) = poisson_solver(ps)%rhsX(ij)
			END IF
			IF ( rY ) THEN
				pen_rhs%Y(i) = poisson_solver(ps)%rhsY(ij)
			END IF
			IF ( rZ ) THEN
				pen_rhs%Z(i) = poisson_solver(ps)%rhsZ(ij)
			END IF

			IF ( lX ) THEN
				pen_lhs%X(i) = poisson_solver(ps)%lhsX(ij)
			END IF
			IF ( lY ) THEN
				pen_lhs%Y(i) = poisson_solver(ps)%lhsY(ij)
			END IF
			IF ( lZ ) THEN
				pen_lhs%Z(i) = poisson_solver(ps)%lhsZ(ij)
			END IF
		END DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Perform ifft
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		CALL pencil_ifft( pen_rhs )
		CALL pencil_ifft( pen_lhs )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store fft field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		DO i = 0,ncell(1)-1
			ij = jn + i
			IF ( rX ) THEN
				poisson_solver(ps)%rhsX(ij) = pen_rhs%X(i)/REAL(ncell(1),MK)
			END IF
			IF ( rY ) THEN
				poisson_solver(ps)%rhsY(ij) = pen_rhs%Y(i)/REAL(ncell(1),MK)
			END IF
			IF ( rZ ) THEN
				poisson_solver(ps)%rhsZ(ij) = pen_rhs%Z(i)/REAL(ncell(1),MK)
			END IF

			IF ( lX ) THEN
				poisson_solver(ps)%lhsX(ij) = pen_lhs%X(i)/REAL(ncell(1),MK)
			END IF
			IF ( lY ) THEN
				poisson_solver(ps)%lhsY(ij) = pen_lhs%Y(i)/REAL(ncell(1),MK)
			END IF
			IF ( lZ ) THEN
				poisson_solver(ps)%lhsZ(ij) = pen_lhs%Z(i)/REAL(ncell(1),MK)
			END IF
		END DO
	END DO

!------------------------------------------------------------------------------!
! Deallocate
!------------------------------------------------------------------------------!
	CALL pencil_deallocate(pen_rhs)
	CALL pencil_deallocate(pen_lhs)

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE poisson_solver_solve2d
