!------------------------------------------------------------------------------!
!  
!  File:         poisson_solver_map.cpp
!  
!  Description:  Performs the MPI mapping of the input communication object
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_map(ps,comm)

USE poisson_solver_communication

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps
	TYPE(class_communication), INTENT(IN) :: comm

!------------------------------------------------------------------------------!
! Variables
!------------------------------------------------------------------------------!
	INTEGER:: nproc, rank, ierr
	INTEGER :: i,j,k
	INTEGER :: p, q, s, sn, sqn, pqs
	INTEGER :: nx, ny
	INTEGER :: iproc, jproc
	INTEGER :: nsend, nrecv
	INTEGER, DIMENSION(3) :: imin, imax
	INTEGER, DIMENSION(3) :: ncell
	INTEGER :: nmap
	LOGICAL :: pack

	LOGICAL :: rX,rY,rZ,lX,lY,lZ,mG

!------------------------------------------------------------------------------!
! Objects
!------------------------------------------------------------------------------!
	TYPE(class_comm_buffer) buffer_send
	TYPE(class_comm_buffer), DIMENSION(:), POINTER :: buffer_recv

!------------------------------------------------------------------------------!
! Get MPI info
!------------------------------------------------------------------------------!
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

!------------------------------------------------------------------------------!
! Find out which fields to map
!------------------------------------------------------------------------------!
	rX = .FALSE.
	rY = .FALSE.
	rZ = .FALSE.
	lX = .FALSE.
	lY = .FALSE.
	lZ = .FALSE.
	mG = .FALSE.

	nmap = 0
	IF( ALLOCATED(poisson_solver(ps)%rhsX) )THEN
		rX = .TRUE.
		nmap = nmap + 1
	END IF
	IF( ALLOCATED(poisson_solver(ps)%rhsY) )THEN
		rY = .TRUE.
		nmap = nmap + 1
	END IF
	IF( ALLOCATED(poisson_solver(ps)%rhsZ) )THEN
		rZ = .TRUE.
		nmap = nmap + 1
	END IF
	IF( ALLOCATED(poisson_solver(ps)%lhsX) )THEN
		lX = .TRUE.
		nmap = nmap + 1
	END IF
	IF( ALLOCATED(poisson_solver(ps)%lhsY) )THEN
		lY = .TRUE.
		nmap = nmap + 1
	END IF
	IF( ALLOCATED(poisson_solver(ps)%lhsZ) )THEN
		lZ = .TRUE.
		nmap = nmap + 1
	END IF
	IF( ALLOCATED(poisson_solver(ps)%mapG) )THEN
		mG = .TRUE.
		nmap = nmap + 1
	END IF

!if(rank .eq. 0)then
!write(*,*) nmap, rX,rY,rZ,lX,lY,lZ,mG
!end if

!------------------------------------------------------------------------------!
! Allocate and pack communication buffers
!------------------------------------------------------------------------------!
	ALLOCATE( buffer_recv(comm%ncomm) )

	DO i = 1,comm%ncomm
!------------------------------------------------------------------------------!
! Self-communication
!------------------------------------------------------------------------------!
		IF ( comm%info(i)%nway .EQ. 0 ) THEN

			imin(1) = comm%info(i)%i2j_min_send(1)
			imin(2) = comm%info(i)%i2j_min_send(2)
			imin(3) = comm%info(i)%i2j_min_send(3)

			imax(1) = comm%info(i)%i2j_max_send(1)
			imax(2) = comm%info(i)%i2j_max_send(2)
			imax(3) = comm%info(i)%i2j_max_send(3)

			nx = comm%info(i)%i2j_ncell_partition_send(1)
			ny = comm%info(i)%i2j_ncell_partition_send(2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Fill receive buffer directly
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			nsend = 2 * nmap * comm%info(i)%i2j_ncell(1) &
			                 * comm%info(i)%i2j_ncell(2) &
			                 * comm%info(i)%i2j_ncell(3)

			ALLOCATE( buffer_recv(i)%val(0:nsend-1) )
			buffer_recv(i)%val = 0.0_MK

			k = 0
			DO s = imin(3),imax(3)
				sn = s * ny
				DO q = imin(2),imax(2)
					sqn = ( sn + q ) * nx
					DO p = imin(1),imax(1)
						pqs = sqn + p

						IF ( rX ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%rhsX(pqs))
							k = k + 1
							buffer_recv(i)%val(k) = AIMAG(poisson_solver(ps)%rhsX(pqs))
							k = k + 1
						END IF
						IF ( rY ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%rhsY(pqs))
							k = k + 1
							buffer_recv(i)%val(k) = AIMAG(poisson_solver(ps)%rhsY(pqs))
							k = k + 1
						END IF
						IF ( rZ ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%rhsZ(pqs))
							k = k + 1
							buffer_recv(i)%val(k) = AIMAG(poisson_solver(ps)%rhsZ(pqs))
							k = k + 1
						END IF

						IF ( lX ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%lhsX(pqs))
							k = k + 1
							buffer_recv(i)%val(k) = AIMAG(poisson_solver(ps)%lhsX(pqs))
							k = k + 1
						END IF
						IF ( lY ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%lhsY(pqs))
							k = k + 1
							buffer_recv(i)%val(k) = AIMAG(poisson_solver(ps)%lhsY(pqs))
							k = k + 1
						END IF
						IF ( lZ ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%lhsZ(pqs))
							k = k + 1
							buffer_recv(i)%val(k) = AIMAG(poisson_solver(ps)%lhsZ(pqs))
							k = k + 1
						END IF

						IF ( mG ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%mapG(pqs))
							k = k + 1
							buffer_recv(i)%val(k) = AIMAG(poisson_solver(ps)%mapG(pqs))
							k = k + 1
						END IF

					END DO
				END DO
			END DO

!------------------------------------------------------------------------------!
! Inter-communication
!------------------------------------------------------------------------------!
		ELSE

			iproc = comm%info(i)%iproc
			jproc = comm%info(i)%jproc

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			pack = .TRUE.
			IF ( rank .EQ. iproc .AND. comm%info(i)%nway .NE. 2 ) THEN

				imin(1) = comm%info(i)%i2j_min_send(1)
				imin(2) = comm%info(i)%i2j_min_send(2)
				imin(3) = comm%info(i)%i2j_min_send(3)

				imax(1) = comm%info(i)%i2j_max_send(1)
				imax(2) = comm%info(i)%i2j_max_send(2)
				imax(3) = comm%info(i)%i2j_max_send(3)

				nsend = 2 * nmap * comm%info(i)%i2j_ncell(1) &
				                 * comm%info(i)%i2j_ncell(2) &
				                 * comm%info(i)%i2j_ncell(3)

				nx = comm%info(i)%i2j_ncell_partition_send(1)
				ny = comm%info(i)%i2j_ncell_partition_send(2)

			ELSE IF ( rank .EQ. jproc .AND. comm%info(i)%nway .GT. 1 ) THEN

				imin(1) = comm%info(i)%j2i_min_send(1)
				imin(2) = comm%info(i)%j2i_min_send(2)
				imin(3) = comm%info(i)%j2i_min_send(3)

				imax(1) = comm%info(i)%j2i_max_send(1)
				imax(2) = comm%info(i)%j2i_max_send(2)
				imax(3) = comm%info(i)%j2i_max_send(3)

				nsend = 2 * nmap * comm%info(i)%j2i_ncell(1) &
				                 * comm%info(i)%j2i_ncell(2) &
				                 * comm%info(i)%j2i_ncell(3)

				nx = comm%info(i)%j2i_ncell_partition_send(1)
				ny = comm%info(i)%j2i_ncell_partition_send(2)
			ELSE
				pack = .FALSE.
			END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Pack send buffer
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			IF (pack) THEN
				ALLOCATE( buffer_send%val(0:nsend-1) )
				buffer_send%val = 0.0_MK
				k = 0

				DO s = imin(3),imax(3)
					sn = s * ny
					DO q = imin(2),imax(2)
						sqn = ( sn + q ) * nx
						DO p = imin(1),imax(1)
							pqs = sqn + p

							IF ( rX ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%rhsX(pqs))
								k = k + 1
								buffer_send%val(k) = AIMAG(poisson_solver(ps)%rhsX(pqs))
								k = k + 1
							END IF
							IF ( rY ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%rhsY(pqs))
								k = k + 1
								buffer_send%val(k) = AIMAG(poisson_solver(ps)%rhsY(pqs))
								k = k + 1
							END IF
							IF ( rZ ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%rhsZ(pqs))
								k = k + 1
								buffer_send%val(k) = AIMAG(poisson_solver(ps)%rhsZ(pqs))
								k = k + 1
							END IF
							IF ( lX ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%lhsX(pqs))
								k = k + 1
								buffer_send%val(k) = AIMAG(poisson_solver(ps)%lhsX(pqs))
								k = k + 1
							END IF
							IF ( lY ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%lhsY(pqs))
								k = k + 1
								buffer_send%val(k) = AIMAG(poisson_solver(ps)%lhsY(pqs))
								k = k + 1
							END IF
							IF ( lZ ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%lhsZ(pqs))
								k = k + 1
								buffer_send%val(k) = AIMAG(poisson_solver(ps)%lhsZ(pqs))
								k = k + 1
							END IF

							IF ( mG ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%mapG(pqs))
								k = k + 1
								buffer_send%val(k) = AIMAG(poisson_solver(ps)%mapG(pqs))
								k = k + 1
							END IF

						END DO
					END DO
				END DO

			END IF

!------------------------------------------------------------------------------!
! Send/Recieve
!------------------------------------------------------------------------------!
			IF ( comm%info(i)%nway .EQ. 1 ) THEN ! only i to j
				IF ( rank .EQ. iproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, jproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. jproc ) THEN ! Receiver
					nrecv = 2 * nmap * comm%info(i)%i2j_ncell(1) & 
					                 * comm%info(i)%i2j_ncell(2) &
					                 * comm%info(i)%i2j_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, iproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( comm%info(i)%nway .EQ. 2 ) THEN ! only j to i
				IF ( rank .EQ. jproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, iproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. iproc ) THEN ! Receiver
					nrecv = 2 * nmap * comm%info(i)%j2i_ncell(1) &
					                 * comm%info(i)%j2i_ncell(2) &
					                 * comm%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( comm%info(i)%nway .EQ. 3 ) THEN ! i to j and j to i
				IF ( rank .EQ. iproc ) THEN
					nrecv = 2 * nmap * comm%info(i)%j2i_ncell(1) &
					                 * comm%info(i)%j2i_ncell(2) &
					                 * comm%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_SENDRECV( buffer_send%val, nsend, MPI_DOUBLE, jproc, iproc, &
					              buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, jproc, &
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				ELSE IF ( rank .EQ. jproc ) THEN
					nrecv = 2 * nmap * comm%info(i)%i2j_ncell(1) &
					                 * comm%info(i)%i2j_ncell(2) &
					                 * comm%info(i)%i2j_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_SENDRECV( buffer_send%val, nsend, MPI_DOUBLE, iproc, jproc, &
					              buffer_recv(i)%val, nrecv, MPI_DOUBLE, iproc, iproc, & 
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				END IF
			END IF ! nway

		END IF

!------------------------------------------------------------------------------!
! De-allocate send buffers
!------------------------------------------------------------------------------!
		IF ( ASSOCIATED(buffer_send%val) ) THEN
			DEALLOCATE( buffer_send%val )
		END IF

!		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	END DO ! ncomm

!------------------------------------------------------------------------------!
! Resize arrays
!------------------------------------------------------------------------------!
	nrecv = comm%partition_recv(rank)%ncell(1) &
	      * comm%partition_recv(rank)%ncell(2) &
	      * comm%partition_recv(rank)%ncell(3)

	IF ( rX ) THEN
		IF( ALLOCATED(poisson_solver(ps)%rhsX) )THEN
			DEALLOCATE(poisson_solver(ps)%rhsX)
		END IF
		ALLOCATE( poisson_solver(ps)%rhsX(0:nrecv-1) )
		poisson_solver(ps)%rhsX = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( rY ) THEN
		IF( ALLOCATED(poisson_solver(ps)%rhsY) )THEN
			DEALLOCATE(poisson_solver(ps)%rhsY)
		END IF
		ALLOCATE( poisson_solver(ps)%rhsY(0:nrecv-1) )
		poisson_solver(ps)%rhsY = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( rZ ) THEN
		IF( ALLOCATED(poisson_solver(ps)%rhsZ) )THEN
			DEALLOCATE(poisson_solver(ps)%rhsZ)
		END IF
		ALLOCATE( poisson_solver(ps)%rhsZ(0:nrecv-1) )
		poisson_solver(ps)%rhsZ = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF

	IF ( lX ) THEN
		IF( ALLOCATED(poisson_solver(ps)%lhsX) )THEN
			DEALLOCATE(poisson_solver(ps)%lhsX)
		END IF
		ALLOCATE( poisson_solver(ps)%lhsX(0:nrecv-1) )
		poisson_solver(ps)%lhsX = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( lY ) THEN
		IF( ALLOCATED(poisson_solver(ps)%lhsY) )THEN
			DEALLOCATE(poisson_solver(ps)%lhsY)
		END IF
		ALLOCATE( poisson_solver(ps)%lhsY(0:nrecv-1) )
		poisson_solver(ps)%lhsY = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( lZ ) THEN
		IF( ALLOCATED(poisson_solver(ps)%lhsZ) )THEN
			DEALLOCATE(poisson_solver(ps)%lhsZ)
		END IF
		ALLOCATE( poisson_solver(ps)%lhsZ(0:nrecv-1) )
		poisson_solver(ps)%lhsZ = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF
	IF ( mG ) THEN
		IF( ALLOCATED(poisson_solver(ps)%mapG) )THEN
			DEALLOCATE(poisson_solver(ps)%mapG)
		END IF
		ALLOCATE( poisson_solver(ps)%mapG(0:nrecv-1) )
		poisson_solver(ps)%mapG = CMPLX(0.0_MK,0.0_MK,MKC)
	END IF

!------------------------------------------------------------------------------!
! Unpack communication buffer
!------------------------------------------------------------------------------!
	DO i = 1,comm%ncomm

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		pack = .TRUE.
		iproc = comm%info(i)%iproc
		jproc = comm%info(i)%jproc

		IF ( rank .EQ. jproc .AND. comm%info(i)%nway .NE. 2 ) THEN
			imin(1) = comm%info(i)%i2j_min_recv(1)
			imin(2) = comm%info(i)%i2j_min_recv(2)
			imin(3) = comm%info(i)%i2j_min_recv(3)

			imax(1) = comm%info(i)%i2j_max_recv(1)
			imax(2) = comm%info(i)%i2j_max_recv(2)
			imax(3) = comm%info(i)%i2j_max_recv(3)

			nx = comm%info(i)%i2j_ncell_partition_recv(1)
			ny = comm%info(i)%i2j_ncell_partition_recv(2)
		ELSE IF ( rank .EQ. iproc .AND. comm%info(i)%nway .GT. 1 ) THEN
			imin(1) = comm%info(i)%j2i_min_recv(1)
			imin(2) = comm%info(i)%j2i_min_recv(2)
			imin(3) = comm%info(i)%j2i_min_recv(3)

			imax(1) = comm%info(i)%j2i_max_recv(1)
			imax(2) = comm%info(i)%j2i_max_recv(2)
			imax(3) = comm%info(i)%j2i_max_recv(3)

			nx = comm%info(i)%j2i_ncell_partition_recv(1)
			ny = comm%info(i)%j2i_ncell_partition_recv(2)
		ELSE
			pack = .FALSE.
		END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Unpack buffers
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		IF (pack) THEN
			k = 0

			DO s = imin(3),imax(3)
				sn = s * ny
				DO q = imin(2),imax(2)
					sqn = ( sn + q ) * nx
					DO p = imin(1),imax(1)
						pqs = sqn + p

						IF ( rX ) THEN
							poisson_solver(ps)%rhsX(pqs) = CMPLX( buffer_recv(i)%val(k), &
							                                  buffer_recv(i)%val(k+1) ,MKC)
							k = k + 2
						END IF
						IF ( rY ) THEN
							poisson_solver(ps)%rhsY(pqs) = CMPLX( buffer_recv(i)%val(k), &
							                                  buffer_recv(i)%val(k+1) ,MKC)
							k = k + 2
						END IF
						IF ( rZ ) THEN
							poisson_solver(ps)%rhsZ(pqs) = CMPLX( buffer_recv(i)%val(k), &
							                                  buffer_recv(i)%val(k+1) ,MKC)
							k = k + 2
						END IF
						IF ( lX ) THEN
							poisson_solver(ps)%lhsX(pqs) = CMPLX( buffer_recv(i)%val(k), &
							                                  buffer_recv(i)%val(k+1) ,MKC)
							k = k + 2
						END IF
						IF ( lY ) THEN
							poisson_solver(ps)%lhsY(pqs) = CMPLX( buffer_recv(i)%val(k), &
							                                  buffer_recv(i)%val(k+1) ,MKC)
							k = k + 2
						END IF
						IF ( lZ ) THEN
							poisson_solver(ps)%lhsZ(pqs) = CMPLX( buffer_recv(i)%val(k), &
							                                  buffer_recv(i)%val(k+1) ,MKC)
							k = k + 2
						END IF

						IF ( mG ) THEN
							poisson_solver(ps)%mapG(pqs) = CMPLX( buffer_recv(i)%val(k), &
							                                  buffer_recv(i)%val(k+1) ,MKC)
							k = k + 2
						END IF

					END DO
				END DO
			END DO

		END IF

	END DO ! ncomm

!------------------------------------------------------------------------------!
! De-allocate recieve buffers
!------------------------------------------------------------------------------!
	IF ( ASSOCIATED(buffer_recv) ) THEN
		DO i = 1,comm%ncomm
			IF( ASSOCIATED( buffer_recv(i)%val ) ) THEN
				DEALLOCATE( buffer_recv(i)%val )
			END IF
		END DO
		DEALLOCATE( buffer_recv )
	END IF

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE poisson_solver_map
