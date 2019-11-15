!------------------------------------------------------------------------------!
!  
!  File:         pencil.f
!  Description:  Includes the module for pencil
!  
!------------------------------------------------------------------------------!
MODULE poisson_solver_pencil

!------------------------------------------------------------------------------!
! Include modules
!------------------------------------------------------------------------------!

IMPLICIT NONE

#include "precision.h"

!------------------------------------------------------------------------------!
! Define pencil class
!------------------------------------------------------------------------------!
	TYPE class_pencil
		INTEGER :: nfft

		LOGICAL :: bX
		LOGICAL :: bY
		LOGICAL :: bZ

		COMPLEX(MKC), DIMENSION(:), ALLOCATABLE :: X
		COMPLEX(MKC), DIMENSION(:), ALLOCATABLE :: Y
		COMPLEX(MKC), DIMENSION(:), ALLOCATABLE :: Z

		COMPLEX(MKC), DIMENSION(:), ALLOCATABLE :: auxX
		COMPLEX(MKC), DIMENSION(:), ALLOCATABLE :: auxY
		COMPLEX(MKC), DIMENSION(:), ALLOCATABLE :: auxZ
	END TYPE class_pencil

!------------------------------------------------------------------------------!
! Include subroutines
!------------------------------------------------------------------------------!
CONTAINS

#include "pencil/pencil_setup.f"

#include "pencil/pencil_resize.f"

#include "pencil/pencil_fft_shift.f"

#include "pencil/pencil_fft.f"

#include "pencil/pencil_ifft.f"

#include "pencil/pencil_deallocate.f"

END MODULE poisson_solver_pencil
