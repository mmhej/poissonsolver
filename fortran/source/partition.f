!------------------------------------------------------------------------------!
!  
!  File:         partition.f
!  Description:  Includes the module for domain partitions
!  
!------------------------------------------------------------------------------!
MODULE poisson_solver_partition

!------------------------------------------------------------------------------!
! Include modules
!------------------------------------------------------------------------------!

IMPLICIT NONE

#include "precision.h"

!------------------------------------------------------------------------------!
! Define field class
!------------------------------------------------------------------------------!
	TYPE class_partition
		INTEGER, DIMENSION(3)  :: ncell
		INTEGER, DIMENSION(3)  :: icell
		REAL(MK), DIMENSION(3) :: dx
	END TYPE class_partition

!------------------------------------------------------------------------------!
! Include subroutines
!------------------------------------------------------------------------------!
CONTAINS

#include "partition/partition_setup.f"

END MODULE poisson_solver_partition
