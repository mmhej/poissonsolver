!------------------------------------------------------------------------------!
! 
!  File:         main.f
! 
!  Description:  
! 
!------------------------------------------------------------------------------!
PROGRAM periodic2d

USE poisson_solver_module

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Parameters
!------------------------------------------------------------------------------!
	INTEGER, PARAMETER :: MK = KIND(1D0)
	REAL(MK),PARAMETER :: pi = 3.1415926535897932_MK
	REAL(MK),PARAMETER :: c  = 10.0_MK
	REAL(MK),PARAMETER :: r0 = 1.0_MK

	INTEGER,DIMENSION(2),PARAMETER :: domain_ncell = (/ 128, 128 /)
	INTEGER,DIMENSION(2),PARAMETER :: domain_bounds = (/ 1, 1 /)

	REAL(MK),DIMENSION(2),PARAMETER :: domain_xmin = (/ -1.0_MK, -1.0_MK /)
	REAL(MK),DIMENSION(2),PARAMETER :: domain_xmax = (/  1.0_MK,  1.0_MK /)


!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER :: rank, nproc, ierr
	INTEGER :: d, n, i, j, jn, ij, ji, pq
	REAL(MK) :: x, y, r

	INTEGER,DIMENSION(2) :: ncell, icell
	REAL(MK),DIMENSION(2) :: xmin, dx

	REAL(MK),DIMENSION(:),ALLOCATABLE :: Ax,Ay,Az,Bx,By,Bz

	REAL(MK) :: err, error
	REAL(MK) :: nrm, norm
	REAL(MK) :: solX, solY, solZ

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
		WRITE(*,*)"TESTING PERIODIC DOMAIN IN 2D... "
		WRITE(*,*)""
		WRITE(*,*)"Number of processors intiated: ", nproc
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


!------------------------------------------------------------------------------!
! Setup domain
!------------------------------------------------------------------------------!
	dx(1) = ( domain_xmax(1) - domain_xmin(1) )/REAL(domain_ncell(1),MK)
	dx(2) = ( domain_xmax(2) - domain_xmin(2) )/REAL(domain_ncell(2),MK)

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

	CALL poisson_solver_initialise( 1 )
	CALL poisson_solver_setup2d( 1, domain_ncell, domain_bounds, dx, 0 )

!------------------------------------------------------------------------------!
! Get mesh info
!------------------------------------------------------------------------------!
	ncell(1) = poisson_solver(1)%partition(rank)%ncell(1)
	ncell(2) = poisson_solver(1)%partition(rank)%ncell(2)

	icell(1) = poisson_solver(1)%partition(rank)%icell(1)
	icell(2) = poisson_solver(1)%partition(rank)%icell(2)

	xmin(1)  = domain_xmin(1) + dx(1) * REAL(poisson_solver(1)%partition(rank)%icell(1),MK)
	xmin(2)  = domain_xmin(2) + dx(2) * REAL(poisson_solver(1)%partition(rank)%icell(2),MK)

!------------------------------------------------------------------------------!
! Allocate fields
!------------------------------------------------------------------------------!
	ALLOCATE( Ax(0:ncell(1)*ncell(2)-1) )
	Ax = 0.0_MK
	ALLOCATE( Bx(0:ncell(1)*ncell(2)-1) )
	Bx = 0.0_MK

!------------------------------------------------------------------------------!
! Initiate fields
!------------------------------------------------------------------------------!
	DO j = 0, ncell(2)-1
		jn = j * ncell(1)
		y  = xmin(2) + (REAL(j,MK) + 0.5_MK)*dx(2)
		DO i = 0, ncell(1)-1
			ij = jn + i
			x  = xmin(1) + (REAL(i,MK) + 0.5_MK)*dx(1)

			Bx(ij) = 8.0_MK * pi * pi * SIN(2.0_MK*pi*x) * sin(2.0_MK*pi*y)

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

	CALL poisson_solver_push(1,Bx,By,Bz)

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

	CALL poisson_solver_solve2d( 1 )

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

	CALL poisson_solver_pull(1,Ax,Ay,Az,Bx,By,Bz)

!------------------------------------------------------------------------------!
! Calculate error integral
!------------------------------------------------------------------------------!
	err   = 0.0_MK
	error = 0.0_MK
	nrm   = 0.0_MK
	norm  = 0.0_MK

	DO j = 0,ncell(2)-1
		jn = j * ncell(1)
		y  = xmin(2) + (REAL(j) + 0.5_MK)*dx(2)
		DO i = 0,ncell(1)-1
			ij = jn + i
			x  = xmin(1) + (REAL(i) + 0.5_MK)*dx(1)

			solX = SIN(2.0_MK*pi*x) * SIN(2.0_MK*pi*y)

			err = err + dx(1)*dx(2)* ( Ax(ij) - solX )**2
			nrm = nrm + dx(1)*dx(2)* solX**2
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
	WRITE(10,'(A,6(I8),A,3(E20.12),A,3(E20.12),A)') "  <ImageData WholeExtent='", icell(1), icell(1) + ncell(1), icell(2), icell(2) + ncell(2), 0, 1, "' Ghostlevel='0' Origin='", domain_xmin(1), domain_xmin(2), 0.0, "' Spacing='", dx(1), dx(2), 0.0, "'>"
	WRITE(10,'(A,6(I8),A)') "    <Piece Extent='", icell(1), icell(1) + ncell(1), icell(2), icell(2) + ncell(2),  0, 1, "'>"
	WRITE(10,'(A)') "      <PointData>"
	WRITE(10,'(A)') "      </PointData>"
	WRITE(10,'(A)') "      <CellData>"
	WRITE(10,'(A)') "        <DataArray type='Float64' Name='B' NumberOfComponents='1'  format='ascii'>"

	DO ij = 0, ncell(1)*ncell(2)-1
		WRITE(10,'(E20.12)') Bx(ij)
	END DO

	WRITE(10,'(A)') "        </DataArray>"
	WRITE(10,'(A)') "        <DataArray type='Float64' Name='A' NumberOfComponents='1'  format='ascii'>"

	DO ij = 0, ncell(1)*ncell(2)-1
		WRITE(10,'(E20.12)') Ax(ij)
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
		WRITE(11,'(A,6(I8),A,3(E20.12),A,3(E20.12),A)') "<PImageData WholeExtent='", 0 , domain_ncell(1), 0 ,domain_ncell(2), 0, 1, "' Ghostlevel='0' Origin='", domain_xmin(1), domain_xmin(2), 0.0, "' Spacing='", dx(1), dx(2), 0.0, "'>"
		WRITE(11,'(A)') "  <PCellData Vectors='output'>"
		WRITE(11,'(A)') "    <PDataArray type='Float64' Name='B' NumberOfComponents='1' format='appended' offset='0'/>"
		WRITE(11,'(A)') "    <PDataArray type='Float64' Name='A' NumberOfComponents='1' format='appended' offset='0'/>"
		WRITE(11,'(A)') "  </PCellData>"
		CLOSE(11)
	END IF

	DO i = 0,nproc-1
		IF ( rank == i ) THEN
			OPEN(12,FILE = TRIM(outfile),POSITION = 'APPEND')
			WRITE(12,'(A,6(I8),3A)') "  <Piece Extent='", icell(1), icell(1) + ncell(1), icell(2), icell(2) + ncell(2), 0, 1, "' Source='", TRIM(filename) , "'/>"
			IF(rank .EQ. nproc - 1) THEN
				WRITE(12,'(A)') "</PImageData>"
				WRITE(12,'(A)') "</VTKFile>"
			END IF
			CLOSE(12)
		END IF
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	END DO

!------------------------------------------------------------------------------!
! Finalise poisson solver
!------------------------------------------------------------------------------!
	CALL poisson_solver_finalise( 1 )

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
END PROGRAM periodic2d
