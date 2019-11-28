!------------------------------------------------------------------------------!
!  
!  File:         poisson_solver_solve3d.f
!  
!  Description:  Solves the Poisson equation in 3D
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_solve3d(ps)

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps
	REAL(MK), PARAMETER  :: pi = 3.1415926535897932_MK

!------------------------------------------------------------------------------!
! Variables
!------------------------------------------------------------------------------!
	INTEGER :: nproc, rank, ierr
	INTEGER :: i, j, k, kn, kjn, ijk
	INTEGER :: nfft
	INTEGER,DIMENSION(3) :: ncell

	COMPLEX(MKC) :: div_psi
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

!------------------------------------------------------------------------------!
! Find out which fields to map
!------------------------------------------------------------------------------!
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
			lZ = .TRUE.
		ELSE
			WRITE(*,'(2A)') " [poisson_solver_solve2d]: ", &
			     "Error rhs fields have not been pushed correctly for grad option." 
			GO TO 9999
		END IF
	ELSE IF (poisson_solver(ps)%lhs_div) THEN
		IF (rX .AND. rY .AND. rZ) THEN
			lX = .TRUE.
		ELSE
			WRITE(*,'(2A)') " [poisson_solver_solve2d]: ", &
			     "Error rhs fields have not been pushed correctly for div option" 
			GO TO 9999
		END IF
	ELSE IF (poisson_solver(ps)%lhs_curl) THEN
		IF (rX .AND. rY .AND. rZ) THEN
			lX = .TRUE.
			lY = .TRUE.
			lZ = .TRUE.
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
	ncell(3) = poisson_solver(ps)%xpen(rank)%ncell(3)

	CALL pencil_resize( pen_rhs, ncell(1) )

	DO k = 0,ncell(3)-1
		kn = k * ncell(2)
		DO j = 0,ncell(2)-1
			kjn = (kn + j) * ncell(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store mesh array in x-pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO i = 0,ncell(1)-1
				ijk = kjn + i
				IF ( rX ) THEN
					pen_rhs%X(i) = poisson_solver(ps)%rhsX(ijk)
				END IF
				IF ( rY ) THEN
					pen_rhs%Y(i) = poisson_solver(ps)%rhsY(ijk)
				END IF
				IF ( rZ ) THEN
					pen_rhs%Z(i) = poisson_solver(ps)%rhsZ(ijk)
				END IF
			END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Fourier transform pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			CALL pencil_fft( pen_rhs )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store pencil in fft array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO i = 0,ncell(1)-1
				ijk = kjn + i
				IF ( rX ) THEN
					poisson_solver(ps)%rhsX(ijk) = pen_rhs%X(i)
				END IF
				IF ( rY ) THEN
					poisson_solver(ps)%rhsY(ijk) = pen_rhs%Y(i)
				END IF
				IF ( rZ ) THEN
					poisson_solver(ps)%rhsZ(ijk) = pen_rhs%Z(i)
				END IF
			END DO
		END DO
	END DO

!------------------------------------------------------------------------------!
! Map to y-pencil partition
!------------------------------------------------------------------------------!
	CALL poisson_solver_map( ps, poisson_solver(ps)%xpen2ypen )

!------------------------------------------------------------------------------!
! FFT y-pencils
!------------------------------------------------------------------------------!
	ncell(1) = poisson_solver(ps)%ypen(rank)%ncell(1)
	ncell(2) = poisson_solver(ps)%ypen(rank)%ncell(2)
	ncell(3) = poisson_solver(ps)%ypen(rank)%ncell(3)

	CALL pencil_resize( pen_rhs, ncell(2) )

	DO k = 0,ncell(3)-1
		kn = k * ncell(2)
		DO i = 0,ncell(1)-1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store mesh array in y-pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO j = 0,ncell(2)-1
				ijk = (kn + j) * ncell(1) + i
				IF ( rX ) THEN
					pen_rhs%X(j) = poisson_solver(ps)%rhsX(ijk)
				END IF
				IF ( rY ) THEN
					pen_rhs%Y(j) = poisson_solver(ps)%rhsY(ijk)
				END IF
				IF ( rZ ) THEN
					pen_rhs%Z(j) = poisson_solver(ps)%rhsZ(ijk)
				END IF
			END DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Fourier transform pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			CALL pencil_fft( pen_rhs )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store pencil in fft array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO j = 0,ncell(2)-1
				ijk = (kn + j) * ncell(1) + i
				IF ( rX ) THEN
					poisson_solver(ps)%rhsX(ijk) = pen_rhs%X(j)
				END IF
				IF ( rY ) THEN
					poisson_solver(ps)%rhsY(ijk) = pen_rhs%Y(j)
				END IF
				IF ( rZ ) THEN
					poisson_solver(ps)%rhsZ(ijk) = pen_rhs%Z(j)
				END IF
			END DO
		END DO
	END DO

!------------------------------------------------------------------------------!
! FFT z-pencils and perform Fourier space operations
!------------------------------------------------------------------------------!
	IF (poisson_solver(ps)%bc(3) .EQ. 0) THEN
		nfft = 2*ncell(3)
	ELSE
		nfft = ncell(3)
	END IF

	CALL pencil_resize( pen_rhs, nfft )
	CALL pencil_resize( pen_lhs, nfft )

	IF ( lX ) THEN
		IF( ALLOCATED(poisson_solver(ps)%lhsX) ) THEN
			IF ( SIZE(poisson_solver(ps)%lhsX) .NE. ncell(1)*ncell(2)*ncell(3) )THEN
				DEALLOCATE( poisson_solver(ps)%lhsX )
				ALLOCATE( poisson_solver(ps)%lhsX( 0:ncell(1)*ncell(2)*ncell(3)-1 ) )
			END IF
		ELSE
			ALLOCATE( poisson_solver(ps)%lhsX( 0:ncell(1)*ncell(2)*ncell(3)-1 ) )
		END IF
		poisson_solver(ps)%lhsX = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( lY ) THEN
		IF( ALLOCATED(poisson_solver(ps)%lhsY) ) THEN
			IF ( SIZE(poisson_solver(ps)%lhsY) .NE. ncell(1)*ncell(2)*ncell(3) )THEN
				DEALLOCATE( poisson_solver(ps)%lhsY )
				ALLOCATE( poisson_solver(ps)%lhsY( 0:ncell(1)*ncell(2)*ncell(3)-1 ) )
			END IF
		ELSE
			ALLOCATE( poisson_solver(ps)%lhsY( 0:ncell(1)*ncell(2)*ncell(3)-1 ) )
		END IF
		poisson_solver(ps)%lhsY = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( lZ ) THEN
		IF( ALLOCATED(poisson_solver(ps)%lhsZ) ) THEN
			IF ( SIZE(poisson_solver(ps)%lhsZ) .NE. ncell(1)*ncell(2)*ncell(3) )THEN
				DEALLOCATE( poisson_solver(ps)%lhsZ )
				ALLOCATE( poisson_solver(ps)%lhsZ( 0:ncell(1)*ncell(2)*ncell(3)-1 ) )
			END IF
		ELSE
			ALLOCATE( poisson_solver(ps)%lhsZ( 0:ncell(1)*ncell(2)*ncell(3)-1 ) )
		END IF
		poisson_solver(ps)%lhsZ = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF

	DO i = 0,ncell(1)-1
		DO j = 0,ncell(2)-1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store mesh array in z-pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO k = 0,ncell(3)-1
				ijk = (k * ncell(2) + j) * ncell(1) + i

				IF ( rX ) THEN
					pen_rhs%X(k) = poisson_solver(ps)%rhsX(ijk)
				END IF
				IF ( rY ) THEN
					pen_rhs%Y(k) = poisson_solver(ps)%rhsY(ijk)
				END IF
				IF ( rZ ) THEN
					pen_rhs%Z(k) = poisson_solver(ps)%rhsZ(ijk)
				END IF
			END DO
! Zero-padding
			IF ( nfft > ncell(3) ) THEN
				DO k = ncell(3),nfft-1
					IF ( rX ) THEN
						pen_rhs%X(k) = CMPLX(0.0_MK,0.0_MK,MKC)
					END IF
					IF ( rY ) THEN
						pen_rhs%Y(k) = CMPLX(0.0_MK,0.0_MK,MKC)
					END IF
					IF ( rZ ) THEN
						pen_rhs%Z(k) = CMPLX(0.0_MK,0.0_MK,MKC)
					END IF
				END DO
			END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! FFT z-pencils
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			CALL pencil_fft( pen_rhs )

!------------------------------------------------------------------------------!
! Perform Fourier space operations
!------------------------------------------------------------------------------!
			DO k = 0,nfft-1
				ijk = (k * ncell(2) + j) * ncell(1) + i

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Regularise rhs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				IF ( poisson_solver(ps)%regularisation .GT. 0 ) THEN
					IF ( rX ) THEN
						pen_rhs%X(k) = pen_rhs%X(k) * poisson_solver(ps)%zeta(ijk)
					END IF
					IF ( rY ) THEN
						pen_rhs%Y(k) = pen_rhs%Y(k) * poisson_solver(ps)%zeta(ijk)
					END IF
					IF ( rZ ) THEN
						pen_rhs%Z(k) = pen_rhs%Z(k) * poisson_solver(ps)%zeta(ijk)
					END IF
				END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Do convolution
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				IF (poisson_solver(ps)%lhs_grad) THEN
					pen_lhs%X(k) = poisson_solver(ps)%ikX(i) * poisson_solver(ps)%G3D(ijk) * pen_rhs%X(k)
					pen_lhs%Y(k) = poisson_solver(ps)%ikY(j) * poisson_solver(ps)%G3D(ijk) * pen_rhs%Y(k)
					pen_lhs%Z(k) = poisson_solver(ps)%ikZ(k) * poisson_solver(ps)%G3D(ijk) * pen_rhs%Z(k)
				ELSE IF (poisson_solver(ps)%lhs_div) THEN
					pen_lhs%X(k) = poisson_solver(ps)%G3D(ijk) * &
					                          ( poisson_solver(ps)%ikX(i) * pen_rhs%X(k) &
					                          + poisson_solver(ps)%ikY(j) * pen_rhs%Y(k) &
					                          + poisson_solver(ps)%ikZ(k) * pen_rhs%Z(k) )
				ELSE IF (poisson_solver(ps)%lhs_curl) THEN
					pen_lhs%X(k) = poisson_solver(ps)%G3D(ijk) * &
					                          ( poisson_solver(ps)%ikY(j) * pen_rhs%Z(k) &
					                          - poisson_solver(ps)%ikZ(k) * pen_rhs%Y(k) )
					pen_lhs%Y(k) = poisson_solver(ps)%G3D(ijk) * &
					                          ( poisson_solver(ps)%ikZ(k) * pen_rhs%X(k) &
					                          - poisson_solver(ps)%ikX(i) * pen_rhs%Z(k) )
					pen_lhs%Z(k) = poisson_solver(ps)%G3D(ijk) * &
					                          ( poisson_solver(ps)%ikX(i) * pen_rhs%Y(k) &
					                          - poisson_solver(ps)%ikY(j) * pen_rhs%X(k) )
				ELSE
					IF ( lX ) THEN
						pen_lhs%X(k) = poisson_solver(ps)%G3D(ijk) * pen_rhs%X(k)
					END IF
					IF ( lY ) THEN
						pen_lhs%Y(k) = poisson_solver(ps)%G3D(ijk) * pen_rhs%Y(k)
					END IF
					IF ( lZ ) THEN
						pen_lhs%Z(k) = poisson_solver(ps)%G3D(ijk) * pen_rhs%Z(k)
					END IF
				END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Reproject rhs field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				IF ( poisson_solver(ps)%rhs_reproject .AND. rX .AND. rY .AND. rZ ) THEN
					div_psi = poisson_solver(ps)%G3D(ijk) * &
					                           ( poisson_solver(ps)%ikX(i) * pen_rhs%X(k) &
					                           + poisson_solver(ps)%ikY(j) * pen_rhs%Y(k) &
					                           + poisson_solver(ps)%ikZ(k) * pen_rhs%Z(k) )
					pen_rhs%X(k) = pen_rhs%X(k) + poisson_solver(ps)%ikX(i) * div_psi
					pen_rhs%Y(k) = pen_rhs%Y(k) + poisson_solver(ps)%ikY(j) * div_psi
					pen_rhs%Z(k) = pen_rhs%Z(k) + poisson_solver(ps)%ikZ(k) * div_psi
				END IF

			END DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! IFFT z-pencils
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			CALL pencil_ifft( pen_rhs )
			CALL pencil_ifft( pen_lhs )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store solution
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO k = 0,ncell(3)-1
				ijk = (k * ncell(2) + j) * ncell(1) + i
				IF ( rX ) THEN
					poisson_solver(ps)%rhsX(ijk) = pen_rhs%X(k)/REAL(nfft,MK)
				END IF
				IF ( rY ) THEN
					poisson_solver(ps)%rhsY(ijk) = pen_rhs%Y(k)/REAL(nfft,MK)
				END IF
				IF ( rZ ) THEN
					poisson_solver(ps)%rhsZ(ijk) = pen_rhs%Z(k)/REAL(nfft,MK)
				END IF

				IF ( lX ) THEN
					poisson_solver(ps)%lhsX(ijk) = pen_lhs%X(k)/REAL(nfft,MK)
				END IF
				IF ( lY ) THEN
					poisson_solver(ps)%lhsY(ijk) = pen_lhs%Y(k)/REAL(nfft,MK)
				END IF
				IF ( lZ ) THEN
					poisson_solver(ps)%lhsZ(ijk) = pen_lhs%Z(k)/REAL(nfft,MK)
				END IF

			END DO
		END DO
	END DO

!------------------------------------------------------------------------------!
! IFFT y-pencils
!------------------------------------------------------------------------------!
	CALL pencil_resize( pen_rhs, ncell(2) )
	CALL pencil_resize( pen_lhs, ncell(2) )

	DO k = 0,ncell(3)-1
		kn = k * ncell(2)
		DO i = 0,ncell(1)-1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store mesh array in y-pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO j = 0,ncell(2)-1
				ijk = (kn + j) * ncell(1) + i

				IF ( rX ) THEN
					pen_rhs%X(j) = poisson_solver(ps)%rhsX(ijk)
				END IF
				IF ( rY ) THEN
					pen_rhs%Y(j) = poisson_solver(ps)%rhsY(ijk)
				END IF
				IF ( rZ ) THEN
					pen_rhs%Z(j) = poisson_solver(ps)%rhsZ(ijk)
				END IF

				IF ( lX ) THEN
					pen_lhs%X(j) = poisson_solver(ps)%lhsX(ijk)
				END IF
				IF ( lY ) THEN
					pen_lhs%Y(j) = poisson_solver(ps)%lhsY(ijk)
				END IF
				IF ( lZ ) THEN
					pen_lhs%Z(j) = poisson_solver(ps)%lhsZ(ijk)
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
			DO j = 0,ncell(2)-1
				ijk = (kn + j) * ncell(1) + i

				IF ( rX ) THEN
					poisson_solver(ps)%rhsX(ijk) = pen_rhs%X(j)/REAL(ncell(2),MK)
				END IF
				IF ( rY ) THEN
					poisson_solver(ps)%rhsY(ijk) = pen_rhs%Y(j)/REAL(ncell(2),MK)
				END IF
				IF ( rZ ) THEN
					poisson_solver(ps)%rhsZ(ijk) = pen_rhs%Z(j)/REAL(ncell(2),MK)
				END IF

				IF ( lX ) THEN
					poisson_solver(ps)%lhsX(ijk) = pen_lhs%X(j)/REAL(ncell(2),MK)
				END IF
				IF ( lY ) THEN
					poisson_solver(ps)%lhsY(ijk) = pen_lhs%Y(j)/REAL(ncell(2),MK)
				END IF
				IF ( lZ ) THEN
					poisson_solver(ps)%lhsZ(ijk) = pen_lhs%Z(j)/REAL(ncell(2),MK)
				END IF
			END DO

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
	ncell(3) = poisson_solver(ps)%xpen(rank)%ncell(3)

	CALL pencil_resize( pen_rhs, ncell(1) )
	CALL pencil_resize( pen_lhs, ncell(1) )

	DO k = 0,ncell(3)-1
		kn = k * ncell(2)
		DO j = 0,ncell(2)-1
			kjn = (kn + j) * ncell(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store mesh array in x-pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO i = 0,ncell(1)-1
				ijk = kjn + i

				IF ( rX ) THEN
					pen_rhs%X(i) = poisson_solver(ps)%rhsX(ijk)
				END IF
				IF ( rY ) THEN
					pen_rhs%Y(i) = poisson_solver(ps)%rhsY(ijk)
				END IF
				IF ( rZ ) THEN
					pen_rhs%Z(i) = poisson_solver(ps)%rhsZ(ijk)
				END IF

				IF ( lX ) THEN
					pen_lhs%X(i) = poisson_solver(ps)%lhsX(ijk)
				END IF
				IF ( lY ) THEN
					pen_lhs%Y(i) = poisson_solver(ps)%lhsY(ijk)
				END IF
				IF ( lZ ) THEN
					pen_lhs%Z(i) = poisson_solver(ps)%lhsZ(ijk)
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
				ijk = kjn + i

				IF ( rX ) THEN
					poisson_solver(ps)%rhsX(ijk) = pen_rhs%X(i)/REAL(ncell(1),MK)
				END IF
				IF ( rY ) THEN
					poisson_solver(ps)%rhsY(ijk) = pen_rhs%Y(i)/REAL(ncell(1),MK)
				END IF
				IF ( rZ ) THEN
					poisson_solver(ps)%rhsZ(ijk) = pen_rhs%Z(i)/REAL(ncell(1),MK)
				END IF

				IF ( lX ) THEN
					poisson_solver(ps)%lhsX(ijk) = pen_lhs%X(i)/REAL(ncell(1),MK)
				END IF
				IF ( lY ) THEN
					poisson_solver(ps)%lhsY(ijk) = pen_lhs%Y(i)/REAL(ncell(1),MK)
				END IF
				IF ( lZ ) THEN
					poisson_solver(ps)%lhsZ(ijk) = pen_lhs%Z(i)/REAL(ncell(1),MK)
				END IF
			END DO

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
END SUBROUTINE poisson_solver_solve3d
