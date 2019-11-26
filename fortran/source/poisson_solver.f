!------------------------------------------------------------------------------!
!  
!  File:         poisson_solver.f
!  Description:  Includes the source code for poisson_solver
!  
!------------------------------------------------------------------------------!
MODULE poisson_solver_module

!------------------------------------------------------------------------------!
! Include modules
!------------------------------------------------------------------------------!
USE poisson_solver_partition
USE poisson_solver_communication
USE poisson_solver_pencil

IMPLICIT NONE

#include "precision.h"

!------------------------------------------------------------------------------!
! Define field class
!------------------------------------------------------------------------------!
	TYPE class_poisson_solver
		LOGICAL,PRIVATE :: lhs_grad
		LOGICAL,PRIVATE :: lhs_div
		LOGICAL,PRIVATE :: lhs_curl
		LOGICAL,PRIVATE :: rhs_reproject
		INTEGER,PRIVATE :: regularisation

		INTEGER, DIMENSION(3)  :: ncell
		INTEGER, DIMENSION(3)  :: bc
		REAL(MK), DIMENSION(3) :: dx

		TYPE(class_partition), DIMENSION(:),ALLOCATABLE :: partition

		TYPE(class_partition), DIMENSION(:),ALLOCATABLE :: xpen
		TYPE(class_partition), DIMENSION(:),ALLOCATABLE :: ypen
		TYPE(class_partition), DIMENSION(:),ALLOCATABLE :: zpen

		TYPE(class_communication) :: real2xpen
		TYPE(class_communication) :: xpen2real

		TYPE(class_communication) :: xpen2ypen
		TYPE(class_communication) :: ypen2xpen

		TYPE(class_communication) :: xpen2zpen
		TYPE(class_communication) :: zpen2xpen

		TYPE(class_communication) :: ypen2zpen
		TYPE(class_communication) :: zpen2ypen

		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: rhsX
		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: rhsY
		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: rhsZ

		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: lhsX
		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: lhsY
		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: lhsZ

		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: ikX
		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: ikY
		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: ikZ

		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: mapG
		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: G2D
		COMPLEX(MKC), DIMENSION(:),ALLOCATABLE :: G3D

		REAL(MK), DIMENSION(:),ALLOCATABLE :: zeta

	END TYPE class_poisson_solver

!------------------------------------------------------------------------------!
! Module variables
!------------------------------------------------------------------------------!
	TYPE(class_poisson_solver) :: poisson_solver

!------------------------------------------------------------------------------!
! Module interfaces
!------------------------------------------------------------------------------!
	INTERFACE poisson_solver_push
		MODULE PROCEDURE poisson_solver_push_1pointer3d
		MODULE PROCEDURE poisson_solver_push_3array
	END INTERFACE

	INTERFACE poisson_solver_pull
		MODULE PROCEDURE poisson_solver_pull_1pointer3d
		MODULE PROCEDURE poisson_solver_pull_3array
	END INTERFACE

!------------------------------------------------------------------------------!
! Include subroutines
!------------------------------------------------------------------------------!
CONTAINS

#include "poisson_solver/poisson_solver_setup2d.f"
#include "poisson_solver/poisson_solver_setup3d.f"

#include "poisson_solver/poisson_solver_settings.f"

#include "poisson_solver/poisson_solver_map.f"

#include "poisson_solver/poisson_solver_special_functions.f"

#include "poisson_solver/poisson_solver_greens2d.f"
#include "poisson_solver/poisson_solver_greens3d.f"

#include "poisson_solver/poisson_solver_solve2d.f"
#include "poisson_solver/poisson_solver_solve3d.f"

#include "poisson_solver/poisson_solver_push.f"
#include "poisson_solver/poisson_solver_pull.f"

#include "poisson_solver/poisson_solver_smooth3d.f"

#include "poisson_solver/poisson_solver_finalise.f"

!------------------------------------------------------------------------------!
! Include customised subroutines
!------------------------------------------------------------------------------!

END MODULE poisson_solver_module
