!------------------------------------------------------------------------------!
!  
!  File:         greens2d.cpp
!  
!  Description:  Calculates the 2D Greens function and other convolution kernels
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_greens2d(dom_ncell,dom_bc,dom_dx,reg_order)

USE poisson_solver_partition
USE poisson_solver_communication
USE poisson_solver_pencil

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  DIMENSION(3), INTENT(IN):: dom_ncell
	INTEGER,  DIMENSION(3), INTENT(IN):: dom_bc
	REAL(MK), DIMENSION(3), INTENT(IN):: dom_dx
	INTEGER, INTENT(IN):: reg_order

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	REAL(MK), PARAMETER ::      pi = 3.1415926535897932_MK
	REAL(MK), PARAMETER ::   gamma = 0.5772156649015329_MK

!------------------------------------------------------------------------------!
! Variables
!------------------------------------------------------------------------------!
	INTEGER :: nproc, rank, ierr
	INTEGER :: i, j, jn, ij, n
	INTEGER :: nfft
	INTEGER, DIMENSION(2) :: ncell
	INTEGER, DIMENSION(2) :: icell
	REAL(MK), DIMENSION(2) :: dx
	REAL(MK), DIMENSION(2) :: xmin
	INTEGER :: mreg
	REAL(MK) :: x,y,r,rho
	REAL(MK) :: arg, sum
	REAL(MK) :: sigma
	REAL(MK) :: h, dk

	REAL(MK), DIMENSION(:), POINTER :: c
	REAL(MK) :: C1, C2, c_1_2pi

	TYPE(class_partition), DIMENSION(:),ALLOCATABLE :: xpen_ext
	TYPE(class_partition), DIMENSION(:),ALLOCATABLE :: ypen_ext

	TYPE(class_communication) :: comm_ext

	TYPE(class_pencil) :: pen

!------------------------------------------------------------------------------!
! Get MPI info
!------------------------------------------------------------------------------!
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

	c_1_2pi = 0.5_MK/pi

!==============================================================================!
! Create extended partitions and communications
!==============================================================================!
	ncell(1) = dom_ncell(1)
	ncell(2) = dom_ncell(2)

	dx(1) = dom_dx(1)
	dx(2) = dom_dx(2)

	h = MAXVAL(dx)

	CALL partition_setup( 1, dom_ncell, dom_bc, dom_dx, .TRUE., .TRUE., .TRUE., xpen_ext )
	CALL partition_setup( 2, dom_ncell, dom_bc, dom_dx, .TRUE., .TRUE., .TRUE., ypen_ext )

	CALL communication_setup( xpen_ext, ypen_ext, comm_ext )

!==============================================================================!
! Create spectral derivative operator on y-pencil partition
!==============================================================================!
	ncell(1) = ypen_ext(rank)%ncell(1)
	ncell(2) = ypen_ext(rank)%ncell(2)

	icell(1) = ypen_ext(rank)%icell(1)
	icell(2) = ypen_ext(rank)%icell(2)

	nfft = xpen_ext(rank)%ncell(1)

	IF( .NOT.ALLOCATED(poisson_solver%ikX) )THEN
		dk = 1.0_MK/( REAL( nfft ,MK) * dx(1) )
		ALLOCATE( poisson_solver%ikX( 0:ncell(1)-1 ) )
		poisson_solver%ikX = CMPLX(0.0_MK,0.0_MK,MKC)

		DO i = 0,ncell(1)-1
			IF ( 2*(icell(1)+i) .LT. nfft ) THEN
				poisson_solver%ikX(i) = CMPLX( 0.0_MK, &
					                      2.0_MK * pi * REAL(icell(1)+i,MK) * dk ,MKC)
			ELSE
				poisson_solver%ikX(i) = CMPLX( 0.0_MK, &
					                      2.0_MK * pi * REAL(icell(1)+i-nfft,MK) * dk ,MKC)
			END IF
		END DO
	END IF

	IF( .NOT.ALLOCATED(poisson_solver%ikY) )THEN
		dk = 1.0_MK/( REAL(ncell(2),MK) * dx(2) )
		ALLOCATE( poisson_solver%ikY( 0:ncell(2)-1 ) )
		poisson_solver%ikY = CMPLX(0.0_MK,0.0_MK,MKC)

		DO j = 0,ncell(2)-1
			IF ( 2*j .LT. ncell(2) ) THEN
				poisson_solver%ikY(j) = CMPLX( 0.0_MK, &
					                      2.0_MK * pi * REAL(j,MK) * dk ,MKC)
			ELSE
				poisson_solver%ikY(j) = CMPLX( 0.0_MK, &
					                      2.0_MK * pi * REAL(j-ncell(2),MK) * dk ,MKC)
			END IF
		END DO
	END IF

!==============================================================================!
! Create regularisation kernel
!==============================================================================!
	IF( .NOT.ALLOCATED(poisson_solver%zeta) )THEN
		sigma = 2.0_MK*h
		IF ( MOD(reg_order, 2) .NE. 0 ) THEN
			mreg = reg_order + 1
		ELSE
			mreg = reg_order
		END IF

		IF ( mreg .GT. 0 ) THEN

			ALLOCATE( poisson_solver%zeta( 0:ncell(1)*ncell(2)-1 ) )
			ALLOCATE( c(mreg/2) )

			c(1) = 1.0_MK
			DO n = 2,mreg/2
				c(n) = c(n-1) * 1.0_MK/REAL( n-1, MK )
			END DO

			DO j = 0,ncell(2)-1
				jn = j * ncell(1)
				DO i = 0,ncell(1)-1
					ij = jn + i
					arg = 0.5_MK * sigma**2 * ( AIMAG(poisson_solver%ikX(i))**2 &
						                        + AIMAG(poisson_solver%ikY(j))**2 )
					sum = 0.0_MK
					DO n = 1,mreg/2
						sum = sum + c(n) * arg**(n-1)
					END DO
					poisson_solver%zeta(ij) = sum * EXP( -arg )
				END DO
			END DO
			DEALLOCATE(c)
		END IF
	END IF

!==============================================================================!
! Create Greens functions
!==============================================================================!
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Unbounded domain
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
	IF (dom_bc(1) .EQ. 0 .AND. dom_bc(2) .EQ. 0) THEN

!------------------------------------------------------------------------------!
! Construct pencils
!------------------------------------------------------------------------------!
		CALL pencil_setup( pen, .TRUE., .FALSE., .FALSE.)

!------------------------------------------------------------------------------!
! FFT x-pencils
!------------------------------------------------------------------------------!
		ncell(1) = xpen_ext(rank)%ncell(1)
		ncell(2) = xpen_ext(rank)%ncell(2)

		icell(1) = xpen_ext(rank)%icell(1)
		icell(2) = xpen_ext(rank)%icell(2)

		dx(1) = xpen_ext(rank)%dx(1)
		dx(2) = xpen_ext(rank)%dx(2)

		CALL pencil_resize( pen, ncell(1) )

		ALLOCATE( poisson_solver%mapG(0:ncell(1)*ncell(2)-1) )
		poisson_solver%mapG = CMPLX(0.0_MK,0.0_MK,MKC)

		xmin(1) = - ( 0.5_MK * REAL( xpen_ext(rank)%ncell(1) ,MK) )*dx(1)
		xmin(2) = - ( 0.5_MK * REAL( ypen_ext(rank)%ncell(2) ,MK) )*dx(2)

		sigma = h/pi ! Spectral
		C2    = LOG(2.0_MK*sigma) - gamma

		DO j = 0,ncell(2)-1
			jn = j * ncell(1)
			y  = xmin(2) + REAL( icell(2) + j, MK)*dx(2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store fft mesh array in x-pencil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO i = 0,ncell(1)-1
				x = xmin(1) + REAL( icell(1) + i ,MK)*dx(1)

! Greens function
				r = SQRT(x*x + y*y)
				rho = r/sigma
				pen%X(i) = CMPLX(-c_1_2pi * (bessel_int_J0( rho ) + C2) * dx(1)*dx(2), & 
				                 0.0_MK,MKC)
			END DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! FFT y-pencils
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			CALL pencil_fft_shift( pen )
			CALL pencil_fft( pen )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO i = 0,ncell(1)-1
				ij = jn + i
				poisson_solver%mapG(ij) = pen%X(i)
			END DO

		END DO

!------------------------------------------------------------------------------!
! Map to y-pencils
!------------------------------------------------------------------------------!
		CALL poisson_solver_map( comm_ext )

!------------------------------------------------------------------------------!
! FFT y-pencils
!------------------------------------------------------------------------------!
		ncell(1) = ypen_ext(rank)%ncell(1)
		ncell(2) = ypen_ext(rank)%ncell(2)

		ALLOCATE( poisson_solver%G2D( 0:ncell(1)*ncell(2)-1 ) )
		poisson_solver%G2D = CMPLX(0.0_MK,0.0_MK,MKC)

		CALL pencil_resize( pen, ncell(2) )

		DO i = 0,ncell(1)-1

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store fft mesh array in y-pencil (fft-shifted)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			DO j = 0,ncell(2)-1
				ij = j * ncell(1) + i
				pen%X(j) = poisson_solver%mapG(ij)
			END DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! FFT y-pencils
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			CALL pencil_fft_shift( pen )
			CALL pencil_fft( pen )

!------------------------------------------------------------------------------!
! Store Fourier space Greens function
!------------------------------------------------------------------------------!
			DO j = 0,ncell(2)-1
				ij = j * ncell(1) + i
				poisson_solver%G2D(ij) = pen%X(j)
			END DO
		END DO

!------------------------------------------------------------------------------!
! De-allocate mapG
!------------------------------------------------------------------------------!
		DEALLOCATE( poisson_solver%mapG )

!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Periodic domain
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
	ELSE IF (dom_bc(1) .EQ. 1 .AND. dom_bc(2) .EQ. 1) THEN

!------------------------------------------------------------------------------!
! Calculate Greens function direcly in Fourier space on y-pencil partition
!------------------------------------------------------------------------------!
		ncell(1) = ypen_ext(rank)%ncell(1)
		ncell(2) = ypen_ext(rank)%ncell(2)

		icell(1) = ypen_ext(rank)%icell(1)
		icell(2) = ypen_ext(rank)%icell(2)

		dx(1)    = ypen_ext(rank)%dx(1)
		dx(2)    = ypen_ext(rank)%dx(2)

		ALLOCATE( poisson_solver%G2D( 0:ncell(1)*ncell(2)-1 ) )

		DO j = 0,ncell(2)-1
			jn = j * ncell(1)
			DO i = 0,ncell(1)-1
				ij = jn + i
				IF ( icell(1) .EQ. 0 .AND. i .EQ. 0 .AND. j .EQ. 0 ) THEN
					poisson_solver%G2D(ij) = CMPLX(0.0_MK,0.0_MK,MKC)
				ELSE
					poisson_solver%G2D(ij) = &
					          -1.0_MK/( poisson_solver%ikX(i)*poisson_solver%ikX(i) &
					                  + poisson_solver%ikY(j)*poisson_solver%ikY(j) )
				END IF
			END DO
		END DO

!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Periodic-unbounded domain
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
	ELSE IF (dom_bc(1) .EQ. 1 .AND. dom_bc(2) .EQ. 0) THEN

		ncell(1) = ypen_ext(rank)%ncell(1)
		ncell(2) = ypen_ext(rank)%ncell(2)

		icell(1) = ypen_ext(rank)%icell(1)
		icell(2) = ypen_ext(rank)%icell(2)

		dx(1)    = ypen_ext(rank)%dx(1)
		dx(2)    = ypen_ext(rank)%dx(2)

		xmin(1) = - ( 0.5_MK * REAL( xpen_ext(rank)%ncell(1) ,MK) )*dx(1)
		xmin(2) = - ( 0.5_MK * REAL( ypen_ext(rank)%ncell(2) ,MK) )*dx(2)

!------------------------------------------------------------------------------!
! Calculate 1D free-space Greens function and Fourier transform it
!------------------------------------------------------------------------------!
		CALL pencil_setup( pen, .TRUE., .FALSE., .FALSE.)
		CALL pencil_resize( pen, ncell(2) )

		sigma = dx(2)/pi ! Spectral
		C1 = sigma/pi
		DO j = 0,ncell(2)-1
			y = xmin(2) + REAL( j ,MK)*dx(2)

! Spectral
			rho = ABS(y)/sigma
			pen%X(j) = CMPLX( -C1 * ( sine_int( rho )*rho + COS(rho) ) * dx(2), &
			                  0.0_MK,MKC)
		END DO

		CALL pencil_fft_shift( pen )
		CALL pencil_fft( pen )

!------------------------------------------------------------------------------!
! Construct Greens function in Fourier space on y-pencil partition
!------------------------------------------------------------------------------!
		ALLOCATE( poisson_solver%G2D( 0:ncell(1)*ncell(2)-1 ) )
		poisson_solver%G2D = CMPLX(0.0_MK,0.0_MK,MKC)

		DO j = 0,ncell(2)-1
			jn = j * ncell(1)
			DO i = 0,ncell(1)-1
				ij = jn + i

				IF ( icell(1) .EQ. 0 .AND. i .EQ. 0 ) THEN
					poisson_solver%G2D(ij) = pen%X(j)
				ELSE
					poisson_solver%G2D(ij) = &
					            -1.0_MK/( poisson_solver%ikX(i)*poisson_solver%ikX(i) &
					                    + poisson_solver%ikY(j)*poisson_solver%ikY(j) )
				END IF

			END DO
		END DO

!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Unbounded-periodic domain
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
	ELSE IF (dom_bc(1) .EQ. 0 .AND. dom_bc(2) .EQ. 1) THEN

		ncell(1) = ypen_ext(rank)%ncell(1)
		ncell(2) = ypen_ext(rank)%ncell(2)

		icell(1) = ypen_ext(rank)%icell(1)
		icell(2) = ypen_ext(rank)%icell(2)

		dx(1)    = ypen_ext(rank)%dx(1)
		dx(2)    = ypen_ext(rank)%dx(2)

		nfft     = xpen_ext(rank)%ncell(1)
		xmin(1)  = - ( 0.5_MK * REAL( nfft ,MK) )*dx(1)

!------------------------------------------------------------------------------!
! Calculate 1D free-space Greens function and Fourier transform it
!------------------------------------------------------------------------------!
		CALL pencil_setup( pen, .TRUE., .FALSE., .FALSE.)
		CALL pencil_resize( pen, nfft )

		sigma = dx(1)/pi ! Spectral
		C1 = sigma/pi
		DO i = 0,nfft-1
			x = xmin(1) + REAL( i ,MK)*dx(1)

! Spectral
			rho = ABS(x)/sigma
			pen%X(j) = CMPLX( -C1 * ( sine_int( rho )*rho + COS(rho) ) * dx(1), &
			                  0.0_MK,MKC)
		END DO

		CALL pencil_fft_shift( pen )
		CALL pencil_fft( pen )

!------------------------------------------------------------------------------!
! Construct Greens function in Fourier space on y-pencil partition
!------------------------------------------------------------------------------!
		ALLOCATE( poisson_solver%G2D( 0:ncell(1)*ncell(2)-1 ) )

		DO j = 0,ncell(2)-1
			jn = j * ncell(1)
			DO i = 0,ncell(1)-1
				ij = jn + i

				IF ( j == 0 ) THEN
					poisson_solver%G2D(ij) = pen%X(icell(1) + i)
				ELSE
					poisson_solver%G2D(ij) = &
					           -1.0_MK/( poisson_solver%ikX(i)*poisson_solver%ikX(i) &
					                   + poisson_solver%ikY(j)*poisson_solver%ikY(j) )
				END IF

			END DO
		END DO

	ELSE
		WRITE(*,'(2A)') " [poisson_solver.greens2d]: ", &
		                "Boundary condition configuration unknown."
	END IF

!------------------------------------------------------------------------------!
! Deallocate
!------------------------------------------------------------------------------!
	CALL pencil_deallocate(pen)
	DEALLOCATE(xpen_ext, ypen_ext)

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE poisson_solver_greens2d

