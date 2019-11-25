!------------------------------------------------------------------------------!
!  
!  File:         setup2d.f
!  Description:  Sets up the poisson solver object for 2D simulations
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_setup2d(ncell,bound_cond,dx,reg_order)

USE poisson_solver_partition

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  DIMENSION(2), INTENT(IN) :: ncell
	INTEGER,  DIMENSION(2), INTENT(IN) :: bound_cond
	REAL(MK), DIMENSION(2), INTENT(IN) :: dx
	INTEGER, OPTIONAL, INTENT(IN) :: reg_order

!------------------------------------------------------------------------------!
! Variables
!------------------------------------------------------------------------------!
	INTEGER :: nproc, rank, ierr
	INTEGER :: i

!------------------------------------------------------------------------------!
! Get MPI info
!------------------------------------------------------------------------------!
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

!------------------------------------------------------------------------------!
! Initialise toggles
!------------------------------------------------------------------------------!
	poisson_solver%lhs_grad = .FALSE.
	poisson_solver%lhs_div = .FALSE.
	poisson_solver%lhs_curl = .FALSE.
	poisson_solver%rhs_reproject = .FALSE.
	IF ( PRESENT(reg_order) ) THEN
		poisson_solver%regularisation = reg_order
	ELSE
		poisson_solver%regularisation = 0
	END IF

!------------------------------------------------------------------------------!
! Pass domain info to partitions
!------------------------------------------------------------------------------!
	poisson_solver%ncell(1) = ncell(1);
	poisson_solver%ncell(2) = ncell(2);
	poisson_solver%ncell(3) = 1;

	poisson_solver%bc(1) = bound_cond(1);
	poisson_solver%bc(2) = bound_cond(2);
	poisson_solver%bc(3) = -1;

	poisson_solver%dx(1) = dx(1);
	poisson_solver%dx(2) = dx(2);
	poisson_solver%dx(3) = 0.0_MK;

!------------------------------------------------------------------------------!
! Setup real topology
!------------------------------------------------------------------------------!
	IF( .NOT.ALLOCATED(poisson_solver%partition) )THEN
		CALL partition_setup( 1, poisson_solver%ncell, poisson_solver%bc, &
		                         poisson_solver%dx, .FALSE., .FALSE., .FALSE., &
		                         poisson_solver%partition )
	END IF

!------------------------------------------------------------------------------!
! Setup Greens function
!------------------------------------------------------------------------------!
	CALL poisson_solver_greens2d( poisson_solver%ncell, poisson_solver%bc, &
	                              poisson_solver%dx,poisson_solver%regularisation)

!------------------------------------------------------------------------------!
! Setup pencil topologies
!------------------------------------------------------------------------------!
	CALL partition_setup( 1, poisson_solver%ncell, poisson_solver%bc, &
	                         poisson_solver%dx, .TRUE., .FALSE., .FALSE., &
	                         poisson_solver%xpen )
	CALL partition_setup( 2, poisson_solver%ncell, poisson_solver%bc, &
	                         poisson_solver%dx, .TRUE., .FALSE., .FALSE., &
	                         poisson_solver%ypen )

!------------------------------------------------------------------------------!
! Construct communicators between topologies
!------------------------------------------------------------------------------!
	CALL communication_setup( poisson_solver%partition, poisson_solver%xpen, &
                            poisson_solver%real2xpen )

	CALL communication_setup( poisson_solver%xpen, poisson_solver%partition, &
                            poisson_solver%xpen2real )

	CALL communication_setup( poisson_solver%xpen, poisson_solver%ypen, &
                            poisson_solver%xpen2ypen )

	CALL communication_setup( poisson_solver%ypen, poisson_solver%xpen, &
                            poisson_solver%ypen2xpen )

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE poisson_solver_setup2d

