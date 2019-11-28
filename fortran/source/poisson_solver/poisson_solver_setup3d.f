!------------------------------------------------------------------------------!
!  
!  File:         setup3d.f
!  Description:  Sets up the poisson solver object for 3D simulations
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_setup3d(ps,ncell,bound_cond,dx,reg_order,reproj)

USE poisson_solver_partition

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps
	INTEGER,  DIMENSION(3), INTENT(IN):: ncell
	INTEGER,  DIMENSION(3), INTENT(IN):: bound_cond
	REAL(MK), DIMENSION(3), INTENT(IN):: dx
	INTEGER, OPTIONAL, INTENT(IN) :: reg_order
	INTEGER, OPTIONAL, INTENT(IN) :: reproj

!------------------------------------------------------------------------------!
! Variables
!------------------------------------------------------------------------------!
	INTEGER:: nproc, rank, ierr

!------------------------------------------------------------------------------!
! Get MPI info
!------------------------------------------------------------------------------!
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

!------------------------------------------------------------------------------!
! Settings
!------------------------------------------------------------------------------!
	poisson_solver(ps)%lhs_grad = .FALSE.
	poisson_solver(ps)%lhs_div  = .FALSE.
	poisson_solver(ps)%lhs_curl = .FALSE.

	IF ( PRESENT(reg_order) ) THEN
		poisson_solver(ps)%regularisation = reg_order
	ELSE
		poisson_solver(ps)%regularisation = 0
	END IF
	IF ( PRESENT(reproj) ) THEN
		poisson_solver(ps)%regularisation = reproj
	ELSE
		poisson_solver(ps)%rhs_reproject = .FALSE.
	END IF

!------------------------------------------------------------------------------!
! Pass domain info to partitions
!------------------------------------------------------------------------------!
	poisson_solver(ps)%ncell(1) = ncell(1);
	poisson_solver(ps)%ncell(2) = ncell(2);
	poisson_solver(ps)%ncell(3) = ncell(3);

	poisson_solver(ps)%bc(1) = bound_cond(1);
	poisson_solver(ps)%bc(2) = bound_cond(2);
	poisson_solver(ps)%bc(3) = bound_cond(3);

	poisson_solver(ps)%dx(1) = dx(1);
	poisson_solver(ps)%dx(2) = dx(2);
	poisson_solver(ps)%dx(3) = dx(3);

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
	CALL poisson_solver_greens3d( ps,poisson_solver(ps)%ncell, &
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
	                         poisson_solver(ps)%dx, .TRUE., .TRUE., .FALSE., &
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
END SUBROUTINE poisson_solver_setup3d

