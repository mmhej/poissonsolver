!------------------------------------------------------------------------------!
!  
!  File:         poisson_solver_smooth3d.f
!  
!  Description:  Solves the Poisson equation in 3D
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_smooth3d(ps,h)

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps
	REAL(MK), INTENT(IN) :: h

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	REAL(MK), PARAMETER :: pi = 3.1415926535897932_MK

!------------------------------------------------------------------------------!
! Variables
!------------------------------------------------------------------------------!
	INTEGER :: nproc, rank, ierr
	INTEGER :: i, j, k, kn, kjn, ijk
	INTEGER :: nfft
	INTEGER,DIMENSION(3) :: ncell
	REAL(MK) :: sigma, arg
	COMPLEX(MKC)       :: zeta
	LOGICAL            :: rX,rY,rZ

	TYPE(class_pencil) :: pen_rhs

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

!------------------------------------------------------------------------------!
! Construct pencils
!------------------------------------------------------------------------------!
	CALL pencil_setup( pen_rhs, rX, rY, rZ)

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

	sigma = h
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
! Smooth rhs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

				arg = -0.5_MK * sigma**2 * ( AIMAG(poisson_solver(ps)%ikX(i))**2 &
				                           + AIMAG(poisson_solver(ps)%ikY(j))**2 &
				                           + AIMAG(poisson_solver(ps)%ikZ(k))**2 )

				zeta = CMPLX( EXP(arg) ,0.0_MK ,MKC)
				IF ( rX ) THEN
					pen_rhs%X(k) = pen_rhs%X(k) * zeta
				END IF
				IF ( rY ) THEN
					pen_rhs%Y(k) = pen_rhs%Y(k) * zeta
				END IF
				IF ( rZ ) THEN
					pen_rhs%Z(k) = pen_rhs%Z(k) * zeta
				END IF

			END DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! IFFT z-pencils
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			CALL pencil_ifft( pen_rhs )

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

			END DO
		END DO
	END DO

!------------------------------------------------------------------------------!
! IFFT y-pencils
!------------------------------------------------------------------------------!
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
! Perform ifft
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			CALL pencil_ifft( pen_rhs )

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
! Perform ifft
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			CALL pencil_ifft( pen_rhs )

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

			END DO

		END DO
	END DO

!------------------------------------------------------------------------------!
! Deallocate
!------------------------------------------------------------------------------!
	CALL pencil_deallocate(pen_rhs)

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE poisson_solver_smooth3d
