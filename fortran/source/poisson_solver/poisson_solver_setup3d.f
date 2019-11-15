!------------------------------------------------------------------------------!
!  
!  File:         setup3d.f
!  Description:  Sets up the poisson solver object for 3D simulations
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_setup3d(ncell,bound_cond,dx,reg_order,reproj)

USE poisson_solver_partition

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
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
	poisson_solver%lhs_grad = .FALSE.
	poisson_solver%lhs_div  = .FALSE.
	poisson_solver%lhs_curl = .FALSE.

	IF ( PRESENT(reg_order) ) THEN
		poisson_solver%regularisation = reg_order
	ELSE
		poisson_solver%regularisation = 0
	END IF
	IF ( PRESENT(reproj) ) THEN
		poisson_solver%regularisation = reproj
	ELSE
		poisson_solver%rhs_reproject = .FALSE.
	END IF

!------------------------------------------------------------------------------!
! Pass domain info to partitions
!------------------------------------------------------------------------------!
	poisson_solver%ncell(1) = ncell(1);
	poisson_solver%ncell(2) = ncell(2);
	poisson_solver%ncell(3) = ncell(3);

	poisson_solver%bc(1) = bound_cond(1);
	poisson_solver%bc(2) = bound_cond(2);
	poisson_solver%bc(3) = bound_cond(3);

	poisson_solver%dx(1) = dx(1);
	poisson_solver%dx(2) = dx(2);
	poisson_solver%dx(3) = dx(3);

!------------------------------------------------------------------------------!
! Setup real topology
!------------------------------------------------------------------------------!
	CALL partition_setup( 1, poisson_solver%ncell, poisson_solver%bc, &
	                         poisson_solver%dx, .FALSE., .FALSE., .FALSE., &
	                         poisson_solver%partition )

!------------------------------------------------------------------------------!
! Setup Greens function
!------------------------------------------------------------------------------!
	CALL poisson_solver_greens3d( poisson_solver%ncell,poisson_solver%bc, & 
	                              poisson_solver%dx,poisson_solver%regularisation)

!------------------------------------------------------------------------------!
! Setup pencil topologies
!------------------------------------------------------------------------------!
	CALL partition_setup( 1, poisson_solver%ncell, poisson_solver%bc, &
	                         poisson_solver%dx, .TRUE., .FALSE., .FALSE., &
	                         poisson_solver%xpen )
	CALL partition_setup( 2, poisson_solver%ncell, poisson_solver%bc, &
	                         poisson_solver%dx, .TRUE., .TRUE., .FALSE., &
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
END SUBROUTINE poisson_solver_setup3d

