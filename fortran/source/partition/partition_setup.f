!------------------------------------------------------------------------------!
!  
!  File:         setup.f
!  Description:  Sets up the domain partition
!  
!------------------------------------------------------------------------------!
SUBROUTINE partition_setup( pencil_dir,domain_ncell,bound_cond,dx, &
                            extendX,extendY,extendZ , part)

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER, INTENT(IN) :: pencil_dir
	INTEGER,  DIMENSION(3), INTENT(IN):: domain_ncell
	INTEGER,  DIMENSION(3), INTENT(IN):: bound_cond
	REAL(MK), DIMENSION(3), INTENT(IN):: dx
	LOGICAL, INTENT(IN) :: extendX
	LOGICAL, INTENT(IN) :: extendY
	LOGICAL, INTENT(IN) :: extendZ

	TYPE(class_partition), DIMENSION(:),ALLOCATABLE, INTENT(OUT) :: part

!------------------------------------------------------------------------------!
! Variables
!------------------------------------------------------------------------------!
	INTEGER:: nproc, rank, ierr

	INTEGER:: i,j,k
	INTEGER, DIMENSION(3) :: nsub
	INTEGER, DIMENSION(3) :: dom_ncell
	INTEGER, DIMENSION(3) :: ncell
	INTEGER, DIMENSION(3) :: icell

	INTEGER, DIMENSION(3) :: rem_ncell
	INTEGER :: iproc

!------------------------------------------------------------------------------!
! Get MPI info
!------------------------------------------------------------------------------!
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

!------------------------------------------------------------------------------!
! Partition into pencil partition given the input direction
!------------------------------------------------------------------------------!
	IF (pencil_dir .EQ. 1) THEN
		nsub(1) = 1
		nsub(2) = nproc
		nsub(3) = 1
	ELSE IF (pencil_dir .EQ. 2) THEN
		nsub(1) = nproc
		nsub(2) = 1
		nsub(3) = 1
	ELSE
!		std::cerr << " [poisson_solver.partition_setup3d]: Error." << std::endl
!		exit(EXIT_FAILURE)
	END IF

!------------------------------------------------------------------------------!
! Extend domain in the X-direction if this is unbounded BCs (y-direction 
! extension will be handled on pencils)
!------------------------------------------------------------------------------!
	IF (bound_cond(1) .EQ. 0 .AND. extendX) THEN
		dom_ncell(1) = 2*domain_ncell(1)
	ELSE
		dom_ncell(1) = domain_ncell(1)
	END IF

	IF (bound_cond(2) .EQ. 0 .AND. extendY) THEN
		dom_ncell(2) = 2*domain_ncell(2)
	ELSE
		dom_ncell(2) = domain_ncell(2)
	END IF

	IF (bound_cond(3) .EQ. 0 .AND. extendZ) THEN
		dom_ncell(3) = 2*domain_ncell(3)
	ELSE
		dom_ncell(3) = domain_ncell(3)
	END IF

!------------------------------------------------------------------------------!
! Assign subdomains and keep track of the local index of the subdomain on
! the processor
!------------------------------------------------------------------------------!
	ALLOCATE( part( 0:nproc-1 ) )

	iproc = 0

	icell(1) = 0
	DO i = 0, nsub(1)-1

		icell(2) = 0
		DO j = 0, nsub(2)-1

			icell(3) = 0
			DO k = 0, nsub(3)-1

!------------------------------------------------------------------------------!
! Divide the grid points to subdomain (ghost not included)
! If the grid points does not equally divide with the number of processors
! the remaining grid points are distributed equally over the last processors
!------------------------------------------------------------------------------!
				rem_ncell(1) = MOD(dom_ncell(1), nsub(1))
				rem_ncell(2) = MOD(dom_ncell(2), nsub(2))
				rem_ncell(3) = MOD(dom_ncell(3), nsub(3))

				IF (i .GE. nsub(1)-rem_ncell(1)) THEN
					ncell(1) = dom_ncell(1)/nsub(1) + 1
				ELSE
					ncell(1) = dom_ncell(1)/nsub(1)
				END IF

				IF (j .GE. nsub(2)-rem_ncell(2)) THEN
					ncell(2) = dom_ncell(2)/nsub(2) + 1
				ELSE
					ncell(2) = dom_ncell(2)/nsub(2)
				END IF

				IF (k .GE. nsub(3)-rem_ncell(3)) THEN
					ncell(3) = dom_ncell(3)/nsub(3) + 1
				ELSE
					ncell(3) = dom_ncell(3)/nsub(3)
				END IF

!------------------------------------------------------------------------------!
! Store the partition
!------------------------------------------------------------------------------!
				part(iproc)%ncell(1) = ncell(1)
				part(iproc)%ncell(2) = ncell(2)
				part(iproc)%ncell(3) = ncell(3)

				part(iproc)%icell(1) = icell(1)
				part(iproc)%icell(2) = icell(2)
				part(iproc)%icell(3) = icell(3)

				part(iproc)%dx(1)    = dx(1)
				part(iproc)%dx(2)    = dx(2)
				part(iproc)%dx(3)    = dx(3)

!------------------------------------------------------------------------------!
! Advance processor and cell counters
!------------------------------------------------------------------------------!
				iproc = iproc + 1
				icell(3) = icell(3) + ncell(3)
			END DO
			icell(2) = icell(2) + ncell(2)
		END DO
		icell(1) = icell(1) + ncell(1)
	END DO

!------------------------------------------------------------------------------!
! Output partition
!------------------------------------------------------------------------------!
!	IF (rank .EQ. 0) THEN
!		DO i = 0,nproc-1
!			WRITE(*,*) "rank:  ", i
!			WRITE(*,*) "ncell: ", part(i)%ncell(1), part(i)%ncell(2), part(i)%ncell(3)
!			WRITE(*,*) "icell: ", part(i)%icell(1), part(i)%icell(2), part(i)%icell(3)
!			WRITE(*,*) "dx:    ", part(i)%dx(1), part(i)%dx(2), part(i)%dx(3)
!			WRITE(*,*) "=="
!		END DO
!	END IF

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE partition_setup

