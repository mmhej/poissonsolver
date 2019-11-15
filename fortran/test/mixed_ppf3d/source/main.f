!------------------------------------------------------------------------------!
! 
!  File:         main.f
! 
!  Description:  
! 
!------------------------------------------------------------------------------!
PROGRAM mixed_ppf3d

USE poisson_solver_module

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	INTEGER, PARAMETER :: MK = KIND(1D0)
	REAL(MK),PARAMETER :: PI = 3.1415926535897932_MK
	REAL(MK),PARAMETER :: c  = 10.0_MK
	REAL(MK),PARAMETER :: r0 = 1.0_MK

	INTEGER,DIMENSION(3),PARAMETER :: domain_ncell  = (/ 64,64,64 /)
	INTEGER,DIMENSION(3),PARAMETER :: domain_bounds = (/ 1,1,0 /)

	REAL(MK),DIMENSION(3),PARAMETER :: domain_xmin = (/ -1.0_MK,-1.0_MK,-1.0_MK /)
	REAL(MK),DIMENSION(3),PARAMETER :: domain_xmax = (/  1.0_MK, 1.0_MK, 1.0_MK /)

!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER  :: rank, nproc, ierr
	INTEGER  :: d, n, i, j, k, kn, kjn, ijk
	REAL(MK) :: x, y, z, r

	INTEGER,DIMENSION(3) :: ncell, icell
	REAL(MK),DIMENSION(3) :: xmin, dx, xmax

	REAL(MK),DIMENSION(:),ALLOCATABLE :: Ax,Ay,Az,Bx,By,Bz

	REAL(MK) :: err, error
	REAL(MK) :: nrm, norm
	REAL(MK) :: Amag, Bmag
	REAL(MK) :: solX, solY, solZ
	REAL(MK) :: diffX, diffY, diffZ

	CHARACTER(LEN=256) :: outfile, filename

!------------------------------------------------------------------------------!
! Initiate program
!------------------------------------------------------------------------------!
	ierr = 0

!------------------------------------------------------------------------------!
! Initialize the OpenMPI library
!------------------------------------------------------------------------------!
	CALL MPI_INIT(ierr)

!------------------------------------------------------------------------------!
! Get MPI info
!------------------------------------------------------------------------------!
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

	IF (rank .EQ. 0) THEN
		WRITE(*,*)""
		WRITE(*,*)"TESTING MIXED PERIODIC-UNBOUNDED DOMAIN IN 3D... "
		WRITE(*,*)""
		WRITE(*,*)"Number of processors intiated: ", nproc
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


!------------------------------------------------------------------------------!
! Setup domain
!------------------------------------------------------------------------------!
	dx(1) = ( domain_xmax(1) - domain_xmin(1) )/REAL(domain_ncell(1),MK)
	dx(2) = ( domain_xmax(2) - domain_xmin(2) )/REAL(domain_ncell(2),MK)
	dx(3) = ( domain_xmax(3) - domain_xmin(3) )/REAL(domain_ncell(3),MK)

!------------------------------------------------------------------------------!
! Setup poisson_solver
!------------------------------------------------------------------------------!
#ifdef __verb
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	IF (rank .EQ. 0) THEN
		WRITE(*,*)" [poisson_solver]: Setup"
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

	CALL poisson_solver_setup3d( domain_ncell, domain_bounds, dx, 0 )
	CALL poisson_solver_set_return_curl( .TRUE. ) ! specify lhs operator

!------------------------------------------------------------------------------!
! Get mesh info
!------------------------------------------------------------------------------!
	ncell(1) = poisson_solver%partition(rank)%ncell(1)
	ncell(2) = poisson_solver%partition(rank)%ncell(2)
	ncell(3) = poisson_solver%partition(rank)%ncell(3)

	icell(1) = poisson_solver%partition(rank)%icell(1)
	icell(2) = poisson_solver%partition(rank)%icell(2)
	icell(3) = poisson_solver%partition(rank)%icell(3)

	xmin(1)  = domain_xmin(1) + dx(1) * REAL(poisson_solver%partition(rank)%icell(1),MK)
	xmin(2)  = domain_xmin(2) + dx(2) * REAL(poisson_solver%partition(rank)%icell(2),MK)
	xmin(3)  = domain_xmin(3) + dx(3) * REAL(poisson_solver%partition(rank)%icell(3),MK)

!------------------------------------------------------------------------------!
! Allocate fields
!------------------------------------------------------------------------------!
	ALLOCATE( Ax(0:ncell(1)*ncell(2)*ncell(3)-1) )
	Ax = 0.0_MK
	ALLOCATE( Ay(0:ncell(1)*ncell(2)*ncell(3)-1) )
	Ay = 0.0_MK
	ALLOCATE( Az(0:ncell(1)*ncell(2)*ncell(3)-1) )
	Az = 0.0_MK

	ALLOCATE( Bx(0:ncell(1)*ncell(2)*ncell(3)-1) )
	Bx = 0.0_MK
	ALLOCATE( By(0:ncell(1)*ncell(2)*ncell(3)-1) )
	By = 0.0_MK
	ALLOCATE( Bz(0:ncell(1)*ncell(2)*ncell(3)-1) )
	Bz = 0.0_MK

!------------------------------------------------------------------------------!
! Initiate fields
!------------------------------------------------------------------------------!
	DO k = 0, ncell(3)-1
		kn = k * ncell(2);
		z  = xmin(3) + (REAL(k,MK) + 0.5_MK)*dx(3);
		DO j = 0, ncell(2)-1
			kjn = (kn + j) * ncell(1)
			y  = xmin(2) + (REAL(j,MK) + 0.5_MK)*dx(2)
			DO i = 0, ncell(1)-1
				ijk = kjn + i
				x  = xmin(1) + (REAL(i,MK) + 0.5_MK)*dx(1)

				IF ( ABS(z) .LT. 1.0 ) THEN

					Bx(ijk) = ( EXP( c* z**2 /((z - 1.0_MK)*(z + 1.0_MK)) ) &
					          * SIN( pi * x) * SIN( pi * y) &
					          * (  2.0_MK * pi**2 * z**8    &
					            -  8.0_MK * pi**2 * z**6    &
					            + 12.0_MK * pi**2 * z**4    &
					            -  6.0_MK *     c * z**4    &
					            -  8.0_MK * pi**2 * z**2    &
					            -  4.0_MK *  c**2 * z**2    &
					            +  4.0_MK *     c * z**2    &
					            +  2.0_MK * pi**2           &
					            +  2.0_MK *         c )     &
					          )*( (z - 1.0_MK)**4 * (z + 1.0_MK)**4 )**(-1)
					By(ijk) = 0.0_MK
					Bz(ijk) = 0.0_MK

				ELSE
					Bx(ijk) = 0.0_MK
					By(ijk) = 0.0_MK
					Bz(ijk) = 0.0_MK
				ENDIF

			END DO
		END DO
	END DO

!------------------------------------------------------------------------------!
! Solve Poisson equation
!------------------------------------------------------------------------------!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Push to poisson_solver array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
#ifdef __verb
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	IF (rank .EQ. 0) THEN
		WRITE(*,*)" [poisson_solver]: Push"
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

	CALL poisson_solver_push(Bx,By,Bz)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Solve
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
#ifdef __verb
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	IF (rank .EQ. 0) THEN
		WRITE(*,*)" [poisson_solver]: Solve"
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

	CALL poisson_solver_solve3d()

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Map back to ClientArray
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
#ifdef __verb
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	IF (rank .EQ. 0) THEN
		WRITE(*,*)" [poisson_solver]: Pull"
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

	CALL poisson_solver_pull(Bx,By,Bz,Ax,Ay,Az)

!------------------------------------------------------------------------------!
! Calculate error integral
!------------------------------------------------------------------------------!
	err   = 0.0_MK
	error = 0.0_MK
	nrm   = 0.0_MK
	norm  = 0.0_MK

	DO k = 0, ncell(3)-1
		kn = k * ncell(2);
		z  = xmin(3) + (REAL(k,MK) + 0.5_MK)*dx(3);
		DO j = 0, ncell(2)-1
			kjn = (kn + j) * ncell(1)
			y  = xmin(2) + (REAL(j,MK) + 0.5_MK)*dx(2)
			DO i = 0, ncell(1)-1
				ijk = kjn + i
				x  = xmin(1) + (REAL(i,MK) + 0.5_MK)*dx(1)

				IF (ABS(z) .LT. 1.0_MK) THEN
					solX = 0.0_MK
					solY = -2.0_MK * c * z &
					       * EXP( c * z**2/((z - 1.0_MK)*(z + 1.0_MK)) ) &
					       * SIN( pi * x ) * SIN( pi * y ) &
					       * ( (z - 1.0_MK)**2 * (z + 1.0_MK)**2 )**(-1)
					solZ = - pi * EXP( c* z**2/( (z - 1.0_MK)*(z + 1.0_MK)) ) &
					       * SIN( pi * x ) * COS( pi * y )
				ELSE
					solX = 0.0_MK
					solY = 0.0_MK
					solZ = 0.0_MK
				END IF

				err = err + dx(1)*dx(2)*dx(3)* ( Ax(ijk) - solX )**2
				nrm = nrm + dx(1)*dx(2)*dx(3)* solX**2

				err = err + dx(1)*dx(2)*dx(3)* ( Ay(ijk) - solY )**2
				nrm = nrm + dx(1)*dx(2)*dx(3)* solY**2

				err = err + dx(1)*dx(2)*dx(3)* ( Az(ijk) - solZ )**2
				nrm = nrm + dx(1)*dx(2)*dx(3)* solZ**2
			END DO
		END DO
	END DO

	CALL MPI_REDUCE( err, error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( nrm, norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	IF (rank == 0) THEN
		WRITE(*,'(A,2(E17.6))') "Error: ", dx(1) , SQRT( error/norm )
	END IF


!------------------------------------------------------------------------------!
! Output field
!------------------------------------------------------------------------------!
#ifdef __verb
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	IF (rank .EQ. 0) THEN
		WRITE(*,*)" [poisson_solver]: Output fields"
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Individual .vti files
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
	WRITE(outfile,'(A,I2.2,A)') './output/mesh_P',rank,'.vti'
	OPEN(10,FILE = TRIM(outfile))
	WRITE(10,'(A)') '<?xml version="1.0"?>'

	WRITE(10,'(A)') "<VTKFile type='ImageData' version='0.1' byte_order='LittleEndian'>"
	WRITE(10,'(A,6(I8),A,3(E20.12),A,3(E20.12),A)') "  <ImageData WholeExtent='", icell(1), icell(1) + ncell(1), icell(2), icell(2) + ncell(2), icell(3), icell(3) + ncell(3), "' Ghostlevel='0' Origin='", domain_xmin(1), domain_xmin(2), domain_xmin(3), "' Spacing='", dx(1), dx(2), dx(3), "'>"
	WRITE(10,'(A,6(I8),A)') "    <Piece Extent='", icell(1), icell(1) + ncell(1), icell(2), icell(2) + ncell(2), icell(3), icell(3) + ncell(3), "'>"
	WRITE(10,'(A)') "      <PointData>"
	WRITE(10,'(A)') "      </PointData>"
	WRITE(10,'(A)') "      <CellData>"
	WRITE(10,'(A)') "        <DataArray type='Float64' Name='B' NumberOfComponents='3'  format='ascii'>"

	DO ijk = 0, ncell(1)*ncell(2)*ncell(3)-1
		WRITE(10,'(3(E20.12))') Bx(ijk), By(ijk), Bz(ijk)
	END DO

	WRITE(10,'(A)') "        </DataArray>"
	WRITE(10,'(A)') "        <DataArray type='Float64' Name='A' NumberOfComponents='3'  format='ascii'>"

	DO ijk = 0, ncell(1)*ncell(2)*ncell(3)-1
		WRITE(10,'(3(E20.12))') Ax(ijk), Ay(ijk), Az(ijk)
	END DO

	WRITE(10,'(A)') "        </DataArray>"
	WRITE(10,'(A)') "      </CellData>"
	WRITE(10,'(A)') "    </Piece>"
	WRITE(10,'(A)') "  </ImageData>"
	WRITE(10,'(A)') "</VTKFile>"

	CLOSE(10)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Main .pvti file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
	WRITE(outfile,'(A,I2.2,A)') './output/mesh.pvti'
	WRITE(filename,'(A,I2.2,A)') './mesh_P',rank,'.vti'

	IF (rank .EQ. 0) THEN
		OPEN(11,FILE = TRIM(outfile))
		WRITE(11,'(A)') '<?xml version="1.0"?>'
		WRITE(11,'(A)') "<VTKFile type='PImageData' version='0.1' byte_order='LittleEndian'>"
		WRITE(11,'(A,6(I8),A,3(E20.12),A,3(E20.12),A)') "<PImageData WholeExtent='", 0 , domain_ncell(1), 0 ,domain_ncell(2), 0, domain_ncell(3), "' Ghostlevel='0' Origin='", domain_xmin(1), domain_xmin(2), domain_xmin(3), "' Spacing='", dx(1), dx(2), dx(3), "'>"
		WRITE(11,'(A)') "  <PCellData Vectors='output'>"
		WRITE(11,'(A)') "    <PDataArray type='Float64' Name='B' NumberOfComponents='3' format='appended' offset='0'/>"
		WRITE(11,'(A)') "    <PDataArray type='Float64' Name='A' NumberOfComponents='3' format='appended' offset='0'/>"
		WRITE(11,'(A)') "  </PCellData>"
		CLOSE(11)
	END IF

	DO i = 0,nproc-1
		IF ( rank == i ) THEN
			OPEN(12,FILE = TRIM(outfile),POSITION = 'APPEND')
			WRITE(12,'(A,6(I8),3A)') "  <Piece Extent='", icell(1), icell(1) + ncell(1), icell(2), icell(2) + ncell(2), icell(3), icell(3) + ncell(3), "' Source='", TRIM(filename) , "'/>"
			IF(rank .EQ. nproc - 1) THEN
				WRITE(12,'(A)') "</PImageData>"
				WRITE(12,'(A)') "</VTKFile>"
			END IF
			CLOSE(12)
		END IF
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	END DO

!------------------------------------------------------------------------------!
! Finalise MPI
!------------------------------------------------------------------------------!
  CALL MPI_finalize(ierr)
  GOTO 1111

!------------------------------------------------------------------------------!
! Return·
!------------------------------------------------------------------------------!
 9999 CONTINUE
  CALL MPI_ABORT(MPI_COMM_WORLD,ierr)

 1111 CONTINUE
END PROGRAM mixed_ppf3d
