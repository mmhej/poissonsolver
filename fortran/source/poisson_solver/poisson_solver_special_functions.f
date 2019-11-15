!------------------------------------------------------------------------------!
!  
!  File:         special_functions.f
!  
!  Description:  Evaluation of special algebraic functions
!  
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! Exponential integral function
! This function is taken from: 
! W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery:
! 'Numerical recipes in Fortran: The art of scientific computing' 2.ed, 1986
!------------------------------------------------------------------------------!
FUNCTION exp_int( n, x ) RESULT(y)

IMPLICIT NONE
!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: n
	REAL(MK), INTENT(IN) :: x
	REAL(MK)             :: y

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	REAL(MK), PARAMETER :: eul = 0.5772156649015329_MK
	INTEGER,  PARAMETER :: maxit = 1000
	REAL(MK), PARAMETER :: eps = 1.0e-10_MK
	REAL(MK), PARAMETER :: fpmin = 1.0e-30_MK

!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER  i, j, nm1
	REAL(MK) a,b,c,d,del,fact,h,psi

!------------------------------------------------------------------------------!
! 
!------------------------------------------------------------------------------!
	y = 0.0_MK
	nm1  = n - 1

	IF ( (n .LT. 0) .OR. (x .LT. 0.0_MK) .OR. &
	     ((x .EQ. 0.0_MK) .AND. ((n .EQ. 0) .OR. (n .EQ. 1)))) THEN
		WRITE(*,*) " [exp_int]: Bad argument"
	ELSE IF (n .EQ. 0) THEN
		y = EXP(-x)/x
	ELSE IF (x .EQ. 0.0_MK) THEN
		y = 1.0_MK/REAL(nm1,MK)
	ELSE IF (x .GT. 1.0_MK) THEN
		b = x + REAL(n,MK)
		c = 1.0_MK/fpmin
		d = 1.0_MK/b
		h = d
		DO i = 1,maxit
			a = REAL(-i*(nm1 + i),MK)
			b = b + 2.0_MK
			d = 1.0_MK/(a*d + b)
			c = b + a/c
			del = c*d
			h = h*del
			IF ( ABS(del - 1.0_MK) .LT. eps ) THEN
				y = h*EXP(-x)
				GO TO 9999
			END IF
		END DO
		WRITE(*,*) " [exp_int]: Continued fractions did not converge."
	ELSE
		IF (nm1 .NE. 0) THEN
			y = 1.0_MK/REAL(nm1,MK)
		ELSE
			y = -LOG(x)-eul
		END IF
		fact = 1.0_MK

		DO i = 1,maxit
			fact = -fact*x/REAL(i,MK)
			IF (i .NE. nm1) THEN
				del = -fact/REAL(i-nm1,MK)
			ELSE
				psi = -eul
				DO j = 1,nm1
					psi = psi + 1.0_MK/REAL(j,MK)
				END DO
				del = fact*(-LOG(x) + psi)
			END IF
			y = y + del
			IF ( ABS(del) .LT. ABS(y)*eps ) THEN
				GO TO 9999
			END IF
		END DO
		WRITE(*,*) " [exp_int]: Series did not converge"
	END IF

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
9999 CONTINUE

RETURN
END FUNCTION exp_int

!------------------------------------------------------------------------------!
! Approximation by Chebyshev polynomials to the integral Bessel function 
! of first kind and order 0: Ji0 = int( q^(-1) (1 - J0) )
! Luke, Y. L: Mathematical functions and their approximations (1975) Table 9.3
!------------------------------------------------------------------------------!
FUNCTION bessel_int_J0( x ) RESULT(y)

IMPLICIT NONE
!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	REAL(MK), INTENT(IN) :: x
	REAL(MK)             :: y

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	REAL(MK), PARAMETER :: pi = 3.1415926535897932_MK
	REAL(MK), PARAMETER :: gamma = 0.5772156649015329_MK

	REAL(MK), DIMENSION(18), PARAMETER :: a(0:17) = (/ 1.35105091918187636388_MK, &
	                                                   0.83791030734868376979_MK, &
	                                                  -0.35047963978529462711_MK, &
	                                                   0.12777415867753198659_MK, &
	                                                  -0.02981035698255560990_MK, &
	                                                   0.00455219841169387328_MK, &
	                                                  -0.00048408621967185359_MK, &
	                                                   0.00003780202859916883_MK, &
	                                                  -0.00000225886908506771_MK, &
	                                                   0.00000010664609068423_MK, &
	                                                  -0.00000000408005443149_MK, &
	                                                   0.00000000012909996251_MK, &
	                                                  -0.00000000000343577839_MK, &
	                                                   0.00000000000007799552_MK, &
	                                                  -0.00000000000000152842_MK, &
	                                                   0.00000000000000002612_MK, &
	                                                  -0.00000000000000000039_MK, &
	                                                   0.00000000000000000001_MK /)

	COMPLEX(MKC), DIMENSION(39), PARAMETER :: c(0:38) =  &
	                     (/ ( 0.95360150809738558095_MK,-0.13917925930200001236_MK), &
	                        (-0.05860838853872331670_MK,-0.12902065726135067062_MK), &
	                        (-0.01020283575659856676_MK, 0.01103004348109535741_MK), &
	                        ( 0.00196012704043622581_MK, 0.00051817180856880364_MK), &
	                        (-0.00009574977697756219_MK,-0.00030928210173975681_MK), &
	                        (-0.00003570479477043714_MK, 0.00004647098443047525_MK), &
	                        ( 0.00001169677960430223_MK,-0.00000008198845340928_MK), &
	                        (-0.00000164386246452682_MK,-0.00000191888381006925_MK), &
	                        (-0.00000007415845751760_MK, 0.00000057813667761104_MK), &
	                        ( 0.00000011434387527717_MK,-0.00000008448997773317_MK), &
	                        (-0.00000003600903214141_MK,-0.00000000525612161520_MK), &
	                        ( 0.00000000601257386446_MK, 0.00000000763257790924_MK), &
	                        ( 0.00000000019124656215_MK,-0.00000000268643963177_MK), &
	                        (-0.00000000054892028385_MK, 0.00000000054279949860_MK), &
	                        ( 0.00000000022740445656_MK,-0.00000000001744365343_MK), &
	                        (-0.00000000005671490865_MK,-0.00000000003975692920_MK), &
	                        ( 0.00000000000607510983_MK, 0.00000000002069683990_MK), &
	                        ( 0.00000000000252060520_MK,-0.00000000000639623674_MK), &
	                        (-0.00000000000191255246_MK, 0.00000000000116359235_MK), &
	                        ( 0.00000000000074056501_MK, 0.00000000000006759603_MK), &
	                        (-0.00000000000018950214_MK,-0.00000000000016557337_MK), &
	                        ( 0.00000000000002021389_MK, 0.00000000000008425597_MK), &
	                        ( 0.00000000000001103617_MK,-0.00000000000002824474_MK), &
	                        (-0.00000000000000889993_MK, 0.00000000000000607698_MK), &
	                        ( 0.00000000000000388558_MK,-0.00000000000000003171_MK), &
	                        (-0.00000000000000119200_MK,-0.00000000000000077237_MK), &
	                        ( 0.00000000000000021456_MK, 0.00000000000000048022_MK), &
	                        ( 0.00000000000000002915_MK,-0.00000000000000019502_MK), &
	                        (-0.00000000000000004877_MK, 0.00000000000000005671_MK), &
	                        ( 0.00000000000000002737_MK,-0.00000000000000000862_MK), &
	                        (-0.00000000000000001080_MK,-0.00000000000000000269_MK), &
	                        ( 0.00000000000000000308_MK, 0.00000000000000000309_MK), &
	                        (-0.00000000000000000042_MK,-0.00000000000000000167_MK), &
	                        (-0.00000000000000000020_MK, 0.00000000000000000066_MK), &
	                        ( 0.00000000000000000020_MK,-0.00000000000000000019_MK), &
	                        (-0.00000000000000000011_MK, 0.00000000000000000003_MK), &
	                        ( 0.00000000000000000004_MK, 0.00000000000000000001_MK), &
	                        (-0.00000000000000000001_MK,-0.00000000000000000001_MK), &
	                        ( 0.00000000000000000000_MK, 0.00000000000000000001_MK) /)

!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER :: i
	REAL(MK),DIMENSION(3) :: T(0:3)
	REAL(MK) :: x8
	REAL(MK) :: x5sh
	REAL(MK) :: fac
	COMPLEX(MKC) sum

!------------------------------------------------------------------------------!
! 
!------------------------------------------------------------------------------!
	x8 = 0.125_MK*x
	x5sh = 2.0_MK*5.0_MK/x - 1.0_MK

	IF ( x .LT. 8.0_MK ) THEN
		T(0) = 1.0_MK
		T(1) = x8

		y = T(0) * a(0)
		DO i = 2,35
			T(2) = 2.0_MK*x8*T(1) - T(0)

			IF( MOD(i,2) .EQ. 0 ) THEN
				y = y + a(i/2) * T(2)
			END IF
			T(0) = T(1)
			T(1) = T(2)
		END DO
	ELSE
		T(0) = 1.0_MK
		T(1) = x5sh

		sum = c(0)*T(0)
		sum = sum + c(1)*T(1)
		DO i = 2,38
			T(2) = 2.0_MK*x5sh*T(1) - T(0)
			T(0) = T(1)
			T(1) = T(2)

			sum = sum + c(i) * T(2)
		END DO
		fac = COS(x + 0.25_MK * pi) * REAL( sum ) &
		    - SIN(x + 0.25_MK * pi) * AIMAG( sum )
		y = SQRT( 2.0_MK/(pi*x) )/x * fac + gamma + LOG(0.5_MK*x)
	END IF

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
9999 CONTINUE

RETURN
END FUNCTION bessel_int_J0


!------------------------------------------------------------------------------!
! Sine integral function
! This function is taken from: 
! W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery:
! 'Numerical recipes in Fortran: The art of scientific computing' 2.ed, 1986
!------------------------------------------------------------------------------!
FUNCTION sine_int( x ) RESULT(y)

IMPLICIT NONE
!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	REAL(MK), INTENT(IN) :: x
	REAL(MK)             :: y

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	REAL(MK), PARAMETER :: pi    = 3.1415926535897932_MK
	REAL(MK), PARAMETER :: gamma = 0.5772156649015329_MK

	INTEGER :: maxit = 1000
	REAL(MK), PARAMETER :: eps = 1.0e-10_MK
	REAL(MK), PARAMETER :: fpmin = 1.0e-30_MK
	REAL(MK), PARAMETER :: fpmax = 1.0e30_MK
	REAL(MK), PARAMETER :: tmin = 2.0_MK

!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER i,k

	REAL(MK) :: a,err,fact,sign,sum,sumc,sums,t,term
	COMPLEX(MKC) h,g,b,c,d,del
	LOGICAL odd

!------------------------------------------------------------------------------!
!
!------------------------------------------------------------------------------!
	y = 0.0_MK

	t = ABS(x)
	IF (t .EQ. 0.0_MK) THEN
		y = 0.0_MK
		GO TO 9999
	END IF

! Evaluate continued fraction by modified Lentz’s method (§5.2).
	IF (t .GT. tmin) THEN
		b = CMPLX(1.0_MK,t)
		c = fpmax
		d = 1.0_MK/b
		h = d

		i = 2
		err = 1.0
		DO WHILE ( i .LT. maxit .AND. err .GT. eps )
			a = - (i-1)**2
			b = b + 2.0_MK
			d = 1.0_MK/(a*d+b) ! Denominators cannot be zero.
			c = b + a/c
			del = c*d
			h = h*del

			err = ABS( REAL( del-1.0_MK ) ) + ABS( AIMAG( del-1.0_MK ) )
			i = i + 1
		END DO

		IF (i .GE. maxit) THEN
			WRITE(*,*) " [sine_int]: Continued fraction failed"
		END IF

		g = CMPLX(COS(t),-SIN(t))
		h = g * h

		y = 0.5_MK*pi + AIMAG(h)
! Evaluate both series simultaneously.
	ELSE

! Special case: avoid failure of convergence test because of underflow.
		IF (t < SQRT(fpmin)) THEN
			sumc = 0.0_MK
			sums = t
		ELSE
			sum  = 0.0_MK
			sums = 0.0_MK
			sumc = 0.0_MK
			sign = 1.0_MK
			fact = 1.0_MK
			odd  = .TRUE.

			k = 1
			err = 1.0_MK
			DO WHILE ( k .LE. maxit .AND. err .GT. eps )
				fact = fact*t/k
				term = fact/k
				sum = sum + sign*term
				err = term/ABS(sum)
				IF (odd) THEN
					sign = -sign
					sums = sum
					sum  = sumc
				ELSE
					sumc = sum
					sum  = sums
				END IF

				odd = .NOT.odd
				k = k + 1
			END DO
			IF (k .GE. maxit) THEN
				WRITE(*,*) " [sine_int]: MAXIT exceeded"
			END IF
		END IF

		y = sums
	END IF

	IF (x .LT. 0.0_MK) THEN
		y = -y
	END IF

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
9999 CONTINUE

RETURN
END FUNCTION sine_int


