!------------------------------------------------------------------------------!
!  
!  File:         setup2d.f
!  Description:  Sets up the poisson solver object for 2D simulations
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_setup2d(ps,ncell,bound_cond,dx,reg_order)

USE poisson_solver_partition

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps
	INTEGER,  DIMENSION(2), INTENT(IN) :: ncell
	INTEGER,  DIMENSION(2), INTENT(IN) :: bound_cond
	REAL(MK), DIMENSION(2), INTENT(IN) :: dx
	INTEGER,  OPTIONAL, INTENT(IN)     :: reg_order

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
	poisson_solver(ps)%lhs_grad = .FALSE.
	poisson_solver(ps)%lhs_div = .FALSE.
	poisson_solver(ps)%lhs_curl = .FALSE.
	poisson_solver(ps)%rhs_reproject = .FALSE.
	IF ( PRESENT(reg_order) ) THEN
		poisson_solver(ps)%regularisation = reg_order
	ELSE
		poisson_solver(ps)%regularisation = 0
	END IF

!------------------------------------------------------------------------------!
! Pass domain info to partitions
!------------------------------------------------------------------------------!
	poisson_solver(ps)%ncell(1) = ncell(1);
	poisson_solver(ps)%ncell(2) = ncell(2);
	poisson_solver(ps)%ncell(3) = 1;

	poisson_solver(ps)%bc(1) = bound_cond(1);
	poisson_solver(ps)%bc(2) = bound_cond(2);
	poisson_solver(ps)%bc(3) = -1;

	poisson_solver(ps)%dx(1) = dx(1);
	poisson_solver(ps)%dx(2) = dx(2);
	poisson_solver(ps)%dx(3) = 0.0_MK;

!------------------------------------------------------------------------------!
! Setup real topology
!------------------------------------------------------------------------------!
	IF( .NOT.ALLOCATED(poisson_solver(ps)%partition) )THEN
		CALL partition_setup( 1, poisson_solver(ps)%ncell, poisson_solver(ps)%bc, &
		                         poisson_solver(ps)%dx, .FALSE., .FALSE., .FALSE., &
		                         poisson_solver(ps)%partition )
	END IF

!------------------------------------------------------------------------------!
! Setup Greens function
!------------------------------------------------------------------------------!
	CALL poisson_solver_greens2d( ps, poisson_solver(ps)%ncell, &
	                              poisson_solver(ps)%bc, &
	                              poisson_solver(ps)%dx, &
	                              poisson_solver(ps)%regularisation)

!------------------------------------------------------------------------------!
! Setup pencil topologies
!------------------------------------------------------------------------------!
	CALL partition_setup( 1, poisson_solver(ps)%ncell, poisson_solver(ps)%bc, &
	                         poisson_solver(ps)%dx, .TRUE., .FALSE., .FALSE., &
	                         poisson_solver(ps)%xpen )
	CALL partition_setup( 2, poisson_solver(ps)%ncell, poisson_solver(ps)%bc, &
	                         poisson_solver(ps)%dx, .TRUE., .FALSE., .FALSE., &
	                         poisson_solver(ps)%ypen )

!------------------------------------------------------------------------------!
! Construct communicators between topologies
!------------------------------------------------------------------------------!
	CALL communication_setup( poisson_solver(ps)%partition, poisson_solver(ps)%xpen, &
                            poisson_solver(ps)%real2xpen )

	CALL communication_setup( poisson_solver(ps)%xpen, poisson_solver(ps)%partition, &
                            poisson_solver(ps)%xpen2real )

	CALL communication_setup( poisson_solver(ps)%xpen, poisson_solver(ps)%ypen, &
                            poisson_solver(ps)%xpen2ypen )

	CALL communication_setup( poisson_solver(ps)%ypen, poisson_solver(ps)%xpen, &
                            poisson_solver(ps)%ypen2xpen )

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE poisson_solver_setup2d

