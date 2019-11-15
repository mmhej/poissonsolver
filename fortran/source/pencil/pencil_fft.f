!------------------------------------------------------------------------------!
!  
!  File:         fft.cpp
!  Description:  Glassmans general N FFT algorithm see W. E. Ferguson: 
!                A simple derivation of Glassmans general N fast Fourier 
!                transform, Comp. & Math. With Appls. 1982
!  
!------------------------------------------------------------------------------!
SUBROUTINE pencil_fft( pen )

IMPLICIT NONE

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	TYPE(class_pencil), INTENT(INOUT) :: pen

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	REAL(MK), PARAMETER :: PI = 3.1415926535897932_MK

!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER :: i,j
	INTEGER :: a, b, c
	INTEGER :: ia, ib, ic, jc
	LOGICAL :: inu

	REAL(MK) :: angle
	COMPLEX(MKC) :: delta, omega
	COMPLEX(MKC) :: sumX,sumY,sumZ

!------------------------------------------------------------------------------!
! Glassman FFT algorithm
!------------------------------------------------------------------------------!
	a = 1
	b = pen%nfft
	c = 1
	inu = .TRUE.

	10 CONTINUE

	IF (b .GT. 1) THEN
		GO TO 30
	END IF

	IF (inu) THEN
		GO TO 9999
	END IF

	DO i = 0,pen%nfft-1
		IF ( pen%bX ) THEN
			pen%X(i) = pen%auxX(i)
		END IF
		IF ( pen%bY ) THEN
			pen%Y(i) = pen%auxY(i)
		END IF
		IF ( pen%bZ ) THEN
			pen%Z(i) = pen%auxZ(i)
		END IF
	END DO

	GO TO 9999

	30 CONTINUE

	a = c*a

	DO c = 2,b
		IF ( MOD(b,c) .EQ. 0) THEN
			GO TO 50
		END IF
	END DO

	50 CONTINUE

	b = b/c

	angle = 2.0_MK*PI/REAL(a*c,MK)

	delta = CMPLX(COS(ANGLE),-SIN(ANGLE),MKC)

	omega = CMPLX(1.0_MK, 0.0_MK,MKC)

	IF (inu) THEN 

		DO ic = 1,c
			DO ia = 1,a
				DO ib = 1,b

					IF ( pen%bX ) THEN
						sumX = pen%X(ib + b*(c - 1) + b*c*(ia - 1) - 1)
					END IF
					IF ( pen%bY ) THEN
						sumY = pen%Y(ib + b*(c - 1) + b*c*(ia - 1) - 1)
					END IF
					IF ( pen%bZ ) THEN
						sumZ = pen%Z(ib + b*(c - 1) + b*c*(ia - 1) - 1)
					END IF

					DO j = 2,c
						jc = c + 1 - j
						IF ( pen%bX ) THEN
							sumX = pen%X(ib + b*(jc-1) + b*c*(ia-1) - 1) + omega * sumX
						END IF
						IF ( pen%bY ) THEN
							sumY = pen%Y(ib + b*(jc-1) + b*c*(ia-1) - 1) + omega * sumY
						END IF
						IF ( pen%bZ ) THEN
							sumZ = pen%Z(ib + b*(jc-1) + b*c*(ia-1) - 1) + omega * sumZ
						END IF
					END DO

					IF ( pen%bX ) THEN
						pen%auxX(ib + b*(ia - 1) + b*a*(ic - 1) - 1) = sumX
					END IF
					IF ( pen%bY ) THEN
						pen%auxY(ib + b*(ia - 1) + b*a*(ic - 1) - 1) = sumY
					END IF
					IF ( pen%bZ ) THEN
						pen%auxZ(ib + b*(ia - 1) + b*a*(ic - 1) - 1) = sumZ
					END IF

				END DO
				omega = delta * omega
			END DO
		END DO

	ELSE

		DO ic = 1,c
			DO ia = 1,a
				DO ib = 1,b

					IF ( pen%bX ) THEN
						sumX = pen%auxX(ib + b*(c - 1) + b*c*(ia - 1) - 1)
					END IF
					IF ( pen%bY ) THEN
						sumY = pen%auxY(ib + b*(c - 1) + b*c*(ia - 1) - 1)
					END IF
					IF ( pen%bZ ) THEN
						sumZ = pen%auxZ(ib + b*(c - 1) + b*c*(ia - 1) - 1)
					END IF

					DO j = 2,c
						jc = c + 1 - j
						IF ( pen%bX ) THEN
							sumX = pen%auxX(ib + b*(jc-1) + b*c*(ia-1) - 1) + omega * sumX
						END IF
						IF ( pen%bY ) THEN
							sumY = pen%auxY(ib + b*(jc-1) + b*c*(ia-1) - 1) + omega * sumY
						END IF
						IF ( pen%bZ ) THEN
							sumZ = pen%auxZ(ib + b*(jc-1) + b*c*(ia-1) - 1) + omega * sumZ
						END IF
					END DO

					IF ( pen%bX ) THEN
						pen%X(ib + b*(ia - 1) + b*a*(ic - 1) - 1) = sumX
					END IF
					IF ( pen%bY ) THEN
						pen%Y(ib + b*(ia - 1) + b*a*(ic - 1) - 1) = sumY
					END IF
					IF ( pen%bZ ) THEN
						pen%Z(ib + b*(ia - 1) + b*a*(ic - 1) - 1) = sumZ
					END IF

				END DO
				omega = delta * omega
			END DO
		END DO

	END IF

	inu = .NOT.inu

	GO TO 10

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
9999 CONTINUE

RETURN
END SUBROUTINE pencil_fft
