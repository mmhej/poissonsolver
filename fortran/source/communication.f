!------------------------------------------------------------------------------!
!  
!  File:         communication.f
!  Description:  Includes the module for communication setup
!  
!------------------------------------------------------------------------------!
MODULE poisson_solver_communication

!------------------------------------------------------------------------------!
! Include modules
!------------------------------------------------------------------------------!

IMPLICIT NONE

#include "precision.h"

!------------------------------------------------------------------------------!
! Define communication info class
!------------------------------------------------------------------------------!
	TYPE class_cominfo
		INTEGER :: iproc
		INTEGER :: jproc

		INTEGER :: nway

		INTEGER, DIMENSION(3) :: i2j_ncell
		INTEGER, DIMENSION(3) :: i2j_min_send
		INTEGER, DIMENSION(3) :: i2j_max_send
		INTEGER, DIMENSION(3) :: i2j_min_recv
		INTEGER, DIMENSION(3) :: i2j_max_recv

		INTEGER, DIMENSION(3) :: i2j_ncell_partition_send
		INTEGER, DIMENSION(3) :: i2j_ncell_partition_recv

		INTEGER, DIMENSION(3) :: j2i_ncell
		INTEGER, DIMENSION(3) :: j2i_min_send
		INTEGER, DIMENSION(3) :: j2i_max_send
		INTEGER, DIMENSION(3) :: j2i_min_recv
		INTEGER, DIMENSION(3) :: j2i_max_recv

		INTEGER, DIMENSION(3) :: j2i_ncell_partition_send
		INTEGER, DIMENSION(3) :: j2i_ncell_partition_recv
	END TYPE class_cominfo

!------------------------------------------------------------------------------!
! Define partition info class
!------------------------------------------------------------------------------!
	TYPE class_partinfo
		INTEGER, DIMENSION(3) :: ncell
	END TYPE class_partinfo

!------------------------------------------------------------------------------!
! Define communication class
!------------------------------------------------------------------------------!
	TYPE class_communication
		INTEGER :: ncomm
		TYPE(class_partinfo), DIMENSION(:),ALLOCATABLE :: partition_send
		TYPE(class_partinfo), DIMENSION(:),ALLOCATABLE :: partition_recv
		TYPE(class_cominfo),  DIMENSION(:),ALLOCATABLE :: info
	END TYPE class_communication

!------------------------------------------------------------------------------!
! Define communication buffer class
!------------------------------------------------------------------------------!
	TYPE class_comm_buffer
		REAL(MK),DIMENSION(:),POINTER :: val => NULL()
	END TYPE class_comm_buffer

!------------------------------------------------------------------------------!
! Include subroutines
!------------------------------------------------------------------------------!
CONTAINS

#include "communication/communication_setup.f"

END MODULE poisson_solver_communication
