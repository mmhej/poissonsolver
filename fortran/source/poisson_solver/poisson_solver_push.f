!------------------------------------------------------------------------------!
!  
!  File:         poisson_solver_push.f
!  
!  Description:  Pushes the field arrays to the poisson solver object
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_push_1pointer3d( ps, offset, real_rhs )

IMPLICIT NONE

include 'mpif.h'
!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps
	INTEGER,DIMENSION(3),INTENT(IN)                 :: offset
	REAL(MK),DIMENSION(:,:,:,:),POINTER, INTENT(IN) :: real_rhs

!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER :: i,j,k
	INTEGER :: p, q, s, sn, sqn, pqs
	INTEGER :: nx, ny
	INTEGER :: nproc, rank, ierr
	INTEGER :: iproc, jproc
	INTEGER :: nsend, nrecv
	INTEGER :: nmap
	LOGICAL :: pack
	INTEGER,DIMENSION(3) :: imin, imax

	LOGICAL :: rX,rY,rZ,lX,lY,lZ

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
! Check that the pointer is associated
!------------------------------------------------------------------------------!
	IF ( .NOT.ASSOCIATED(real_rhs) ) THEN
		WRITE(*,'(A)')"[poisson_solver_push]: Error pointer not associated."
		GO TO 9999
	END IF

!------------------------------------------------------------------------------!
! Find out which fields to map
!------------------------------------------------------------------------------!
	rX = .FALSE.
	rY = .FALSE.
	rZ = .FALSE.
	lX = .FALSE.
	lY = .FALSE.
	lZ = .FALSE.

	nmap = 0
	IF ( SIZE(real_rhs,1) .EQ. 1 ) THEN
		rX = .TRUE.
		nmap = 1
	END IF
	IF ( SIZE(real_rhs,1) .EQ. 2 ) THEN
		rX = .TRUE.
		rY = .TRUE.
		nmap = 2
	END IF
	IF ( SIZE(real_rhs,1) .EQ. 3 ) THEN
		rX = .TRUE.
		rY = .TRUE.
		rZ = .TRUE.
		nmap = 3
	END IF

!------------------------------------------------------------------------------!
! Allocate and pack communication buffers
!------------------------------------------------------------------------------!
	ALLOCATE( buffer_recv(poisson_solver(ps)%real2xpen%ncomm) )

	DO i = 1,poisson_solver(ps)%real2xpen%ncomm
!------------------------------------------------------------------------------!
! Self-communication
!------------------------------------------------------------------------------!
		IF ( poisson_solver(ps)%real2xpen%info(i)%nway .EQ. 0 ) THEN

			imin(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(1)
			imin(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(2)
			imin(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(3)

			imax(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(1)
			imax(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(2)
			imax(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(3)

			nx = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_send(1)
			ny = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_send(2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Fill receive buffer directly
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			nsend = nmap * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(1) &
			             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(2) &
			             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(3)

			ALLOCATE( buffer_recv(i)%val(0:nsend-1) )

			k = 0
			DO s = imin(3),imax(3)
				sn = s * ny
				DO q = imin(2),imax(2)
					sqn = ( sn + q ) * nx
					DO p = imin(1),imax(1)
						pqs = sqn + p

						IF ( rX ) THEN
							buffer_recv(i)%val(k) = & 
							            & real_rhs(1,p+offset(1),q+offset(2),s+offset(3))
							k = k + 1
						END IF
						IF ( rY ) THEN
							buffer_recv(i)%val(k) = & 
							            & real_rhs(2,p+offset(1),q+offset(2),s+offset(3))
							k = k + 1
						END IF
						IF ( rZ ) THEN
							buffer_recv(i)%val(k) = & 
							            & real_rhs(3,p+offset(1),q+offset(2),s+offset(3))
							k = k + 1
						END IF

					END DO
				END DO
			END DO

!------------------------------------------------------------------------------!
! Inter-communication
!------------------------------------------------------------------------------!
		ELSE

			iproc = poisson_solver(ps)%real2xpen%info(i)%iproc
			jproc = poisson_solver(ps)%real2xpen%info(i)%jproc

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			pack = .TRUE.
			IF ( rank .EQ. iproc .AND. &
			     poisson_solver(ps)%real2xpen%info(i)%nway .NE. 2 ) THEN
				imin(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(1)
				imin(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(2)
				imin(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(3)

				imax(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(1)
				imax(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(2)
				imax(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(3)

				nsend = nmap * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(1) &
				             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(2) &
				             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(3)

				nx = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_send(1)
				ny = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_send(2)
			ELSE IF ( rank .EQ. jproc .AND. &
			          poisson_solver(ps)%real2xpen%info(i)%nway .GT. 1 ) THEN
				imin(1) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_send(1)
				imin(2) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_send(2)
				imin(3) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_send(3)

				imax(1) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_send(1)
				imax(2) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_send(2)
				imax(3) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_send(3)

				nsend = nmap * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(1) &
				             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(2) &
				             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(3)

				nx = poisson_solver(ps)%real2xpen%info(i)%j2i_ncell_partition_send(1)
				ny = poisson_solver(ps)%real2xpen%info(i)%j2i_ncell_partition_send(2)
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
								buffer_send%val(k) = & 
								          & real_rhs(1,p+offset(1),q+offset(2),s+offset(3))
								k = k + 1
							END IF
							IF ( rY ) THEN
								buffer_send%val(k) = & 
								          & real_rhs(2,p+offset(1),q+offset(2),s+offset(3))
								k = k + 1
							END IF
							IF ( rZ ) THEN
								buffer_send%val(k) = & 
								          & real_rhs(3,p+offset(1),q+offset(2),s+offset(3))
								k = k + 1
							END IF

						END DO
					END DO
				END DO

			END IF

!------------------------------------------------------------------------------!
! Send/Recieve
!------------------------------------------------------------------------------!
			IF ( poisson_solver(ps)%real2xpen%info(i)%nway .EQ. 1 ) THEN ! only i to j
				IF ( rank .EQ. iproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, jproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. jproc ) THEN ! Receiver
					nrecv = nmap * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(1) &
					             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(2) &
					             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, iproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( poisson_solver(ps)%real2xpen%info(i)%nway .EQ. 2 ) THEN ! only j to i
				IF ( rank .EQ. jproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, iproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. iproc ) THEN ! Receiver
					nrecv = nmap * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(1) &
					             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(2) &
					             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( poisson_solver(ps)%real2xpen%info(i)%nway .EQ. 3 ) THEN ! i to j and j to i
				IF ( rank .EQ. iproc ) THEN
					nrecv = nmap * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(1) &
					             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(2) &
					             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_SENDRECV( buffer_send%val, nsend, MPI_DOUBLE, jproc, iproc, &
					              buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, jproc, &
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				ELSE IF ( rank .EQ. jproc ) THEN
					nrecv = nmap * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(1) &
					             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(2) &
					             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_SENDRECV( buffer_send%val, nsend, MPI_DOUBLE, iproc, jproc, &
					              buffer_recv(i)%val, nrecv, MPI_DOUBLE, iproc, iproc, & 
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				END IF
			END IF ! nway

		END IF

!		MPI_Barrier(MPI_COMM_WORLD)

!------------------------------------------------------------------------------!
! De-allocate send buffers
!------------------------------------------------------------------------------!
		IF ( ASSOCIATED(buffer_send%val) ) THEN
			DEALLOCATE( buffer_send%val )
		END IF

	END DO ! ncomm

!------------------------------------------------------------------------------!
! Resize arrays
!------------------------------------------------------------------------------!
	nrecv = poisson_solver(ps)%real2xpen%partition_recv(rank)%ncell(1) &
	      * poisson_solver(ps)%real2xpen%partition_recv(rank)%ncell(2) &
	      * poisson_solver(ps)%real2xpen%partition_recv(rank)%ncell(3)

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


!------------------------------------------------------------------------------!
! Unpack communication buffer
!------------------------------------------------------------------------------!
	DO i = 1,poisson_solver(ps)%real2xpen%ncomm

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		pack = .TRUE.
		iproc = poisson_solver(ps)%real2xpen%info(i)%iproc
		jproc = poisson_solver(ps)%real2xpen%info(i)%jproc

		IF ( rank .EQ. jproc .AND. poisson_solver(ps)%real2xpen%info(i)%nway .NE. 2 ) THEN

			imin(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_recv(1)
			imin(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_recv(2)
			imin(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_recv(3)

			imax(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_recv(1)
			imax(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_recv(2)
			imax(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_recv(3)

			nx = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_recv(1)
			ny = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_recv(2)

		ELSE IF ( rank .EQ. iproc .AND. poisson_solver(ps)%real2xpen%info(i)%nway .GT. 1 ) THEN

			imin(1) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_recv(1)
			imin(2) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_recv(2)
			imin(3) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_recv(3)

			imax(1) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_recv(1)
			imax(2) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_recv(2)
			imax(3) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_recv(3)

			nx = poisson_solver(ps)%real2xpen%info(i)%j2i_ncell_partition_recv(1)
			ny = poisson_solver(ps)%real2xpen%info(i)%j2i_ncell_partition_recv(2)
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
							poisson_solver(ps)%rhsX(pqs) = CMPLX( buffer_recv(i)%val(k), 0.0 ,MKC)
							k = k + 1
						END IF
						IF ( rY ) THEN
							poisson_solver(ps)%rhsY(pqs) = CMPLX( buffer_recv(i)%val(k), 0.0 ,MKC)
							k = k + 1
						END IF
						IF ( rZ ) THEN
							poisson_solver(ps)%rhsZ(pqs) = CMPLX( buffer_recv(i)%val(k), 0.0 ,MKC)
							k = k + 1
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
		DO i = 1,poisson_solver(ps)%real2xpen%ncomm
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
END SUBROUTINE poisson_solver_push_1pointer3d







SUBROUTINE poisson_solver_push_3array( ps, real_rhsX, real_rhsY, real_rhsZ )

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps
	REAL(MK),DIMENSION(:),ALLOCATABLE, INTENT(IN) :: real_rhsX
	REAL(MK),DIMENSION(:),ALLOCATABLE, INTENT(IN) :: real_rhsY
	REAL(MK),DIMENSION(:),ALLOCATABLE, INTENT(IN) :: real_rhsZ

!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER :: i,j,k
	INTEGER :: p, q, s, sn, sqn, pqs
	INTEGER :: nx, ny
	INTEGER :: nproc, rank, ierr
	INTEGER :: iproc, jproc
	INTEGER :: nsend, nrecv
	INTEGER :: nmap
	LOGICAL :: pack
	INTEGER,DIMENSION(3) :: imin, imax

	LOGICAL :: rX,rY,rZ,lX,lY,lZ

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

	nmap = 0
	IF ( ALLOCATED(real_rhsX) ) THEN
		rX = .TRUE.
		nmap = nmap + 1
	END IF
	IF ( ALLOCATED(real_rhsY) ) THEN
		rY = .TRUE.
		nmap = nmap + 1
	END IF
	IF ( ALLOCATED(real_rhsZ) ) THEN
		rZ = .TRUE.
		nmap = nmap + 1
	END IF

!------------------------------------------------------------------------------!
! Allocate and pack communication buffers
!------------------------------------------------------------------------------!
	ALLOCATE( buffer_recv(poisson_solver(ps)%real2xpen%ncomm) )

	DO i = 1,poisson_solver(ps)%real2xpen%ncomm
!------------------------------------------------------------------------------!
! Self-communication
!------------------------------------------------------------------------------!
		IF ( poisson_solver(ps)%real2xpen%info(i)%nway .EQ. 0 ) THEN

			imin(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(1)
			imin(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(2)
			imin(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(3)

			imax(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(1)
			imax(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(2)
			imax(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(3)

			nx = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_send(1)
			ny = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_send(2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Fill receive buffer directly
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			nsend = nmap * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(1) &
			             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(2) &
			             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(3)

			ALLOCATE( buffer_recv(i)%val(0:nsend-1) )

			k = 0
			DO s = imin(3),imax(3)
				sn = s * ny
				DO q = imin(2),imax(2)
					sqn = ( sn + q ) * nx
					DO p = imin(1),imax(1)
						pqs = sqn + p

						IF ( rX ) THEN
							buffer_recv(i)%val(k) = real_rhsX(pqs)
							k = k + 1
						END IF
						IF ( rY ) THEN
							buffer_recv(i)%val(k) = real_rhsY(pqs)
							k = k + 1
						END IF
						IF ( rZ ) THEN
							buffer_recv(i)%val(k) = real_rhsZ(pqs)
							k = k + 1
						END IF

					END DO
				END DO
			END DO

!------------------------------------------------------------------------------!
! Inter-communication
!------------------------------------------------------------------------------!
		ELSE

			iproc = poisson_solver(ps)%real2xpen%info(i)%iproc
			jproc = poisson_solver(ps)%real2xpen%info(i)%jproc

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			pack = .TRUE.
			IF ( rank .EQ. iproc .AND. &
			     poisson_solver(ps)%real2xpen%info(i)%nway .NE. 2 ) THEN
				imin(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(1)
				imin(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(2)
				imin(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_send(3)

				imax(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(1)
				imax(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(2)
				imax(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_send(3)

				nsend = nmap * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(1) &
				             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(2) &
				             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(3)

				nx = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_send(1)
				ny = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_send(2)
			ELSE IF ( rank .EQ. jproc .AND. &
			          poisson_solver(ps)%real2xpen%info(i)%nway .GT. 1 ) THEN
				imin(1) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_send(1)
				imin(2) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_send(2)
				imin(3) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_send(3)

				imax(1) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_send(1)
				imax(2) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_send(2)
				imax(3) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_send(3)

				nsend = nmap * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(1) &
				             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(2) &
				             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(3)

				nx = poisson_solver(ps)%real2xpen%info(i)%j2i_ncell_partition_send(1)
				ny = poisson_solver(ps)%real2xpen%info(i)%j2i_ncell_partition_send(2)
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
								buffer_send%val(k) = real_rhsX(pqs)
								k = k + 1
							END IF
							IF ( rY ) THEN
								buffer_send%val(k) = real_rhsY(pqs)
								k = k + 1
							END IF
							IF ( rZ ) THEN
								buffer_send%val(k) = real_rhsZ(pqs)
								k = k + 1
							END IF

						END DO
					END DO
				END DO

			END IF

!------------------------------------------------------------------------------!
! Send/Recieve
!------------------------------------------------------------------------------!
			IF ( poisson_solver(ps)%real2xpen%info(i)%nway .EQ. 1 ) THEN ! only i to j
				IF ( rank .EQ. iproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, jproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. jproc ) THEN ! Receiver
					nrecv = nmap * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(1) &
					             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(2) &
					             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, iproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( poisson_solver(ps)%real2xpen%info(i)%nway .EQ. 2 ) THEN ! only j to i
				IF ( rank .EQ. jproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, iproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. iproc ) THEN ! Receiver
					nrecv = nmap * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(1) &
					             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(2) &
					             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( poisson_solver(ps)%real2xpen%info(i)%nway .EQ. 3 ) THEN ! i to j and j to i
				IF ( rank .EQ. iproc ) THEN
					nrecv = nmap * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(1) &
					             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(2) &
					             * poisson_solver(ps)%real2xpen%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_SENDRECV( buffer_send%val, nsend, MPI_DOUBLE, jproc, iproc, &
					              buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, jproc, &
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				ELSE IF ( rank .EQ. jproc ) THEN
					nrecv = nmap * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(1) &
					             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(2) &
					             * poisson_solver(ps)%real2xpen%info(i)%i2j_ncell(3)
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

!		MPI_Barrier(MPI_COMM_WORLD)
	END DO ! ncomm


!------------------------------------------------------------------------------!
! Resize arrays
!------------------------------------------------------------------------------!
	nrecv = poisson_solver(ps)%real2xpen%partition_recv(rank)%ncell(1) &
	      * poisson_solver(ps)%real2xpen%partition_recv(rank)%ncell(2) &
	      * poisson_solver(ps)%real2xpen%partition_recv(rank)%ncell(3)

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


!------------------------------------------------------------------------------!
! Unpack communication buffer
!------------------------------------------------------------------------------!
	DO i = 1,poisson_solver(ps)%real2xpen%ncomm

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		pack = .TRUE.
		iproc = poisson_solver(ps)%real2xpen%info(i)%iproc
		jproc = poisson_solver(ps)%real2xpen%info(i)%jproc

		IF ( rank .EQ. jproc .AND. poisson_solver(ps)%real2xpen%info(i)%nway .NE. 2 ) THEN

			imin(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_recv(1)
			imin(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_recv(2)
			imin(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_min_recv(3)

			imax(1) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_recv(1)
			imax(2) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_recv(2)
			imax(3) = poisson_solver(ps)%real2xpen%info(i)%i2j_max_recv(3)

			nx = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_recv(1)
			ny = poisson_solver(ps)%real2xpen%info(i)%i2j_ncell_partition_recv(2)

		ELSE IF ( rank .EQ. iproc .AND. poisson_solver(ps)%real2xpen%info(i)%nway .GT. 1 ) THEN

			imin(1) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_recv(1)
			imin(2) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_recv(2)
			imin(3) = poisson_solver(ps)%real2xpen%info(i)%j2i_min_recv(3)

			imax(1) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_recv(1)
			imax(2) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_recv(2)
			imax(3) = poisson_solver(ps)%real2xpen%info(i)%j2i_max_recv(3)

			nx = poisson_solver(ps)%real2xpen%info(i)%j2i_ncell_partition_recv(1)
			ny = poisson_solver(ps)%real2xpen%info(i)%j2i_ncell_partition_recv(2)
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
							poisson_solver(ps)%rhsX(pqs) = CMPLX( buffer_recv(i)%val(k), 0.0 ,MKC)
							k = k + 1
						END IF
						IF ( rY ) THEN
							poisson_solver(ps)%rhsY(pqs) = CMPLX( buffer_recv(i)%val(k), 0.0 ,MKC)
							k = k + 1
						END IF
						IF ( rZ ) THEN
							poisson_solver(ps)%rhsZ(pqs) = CMPLX( buffer_recv(i)%val(k), 0.0 ,MKC)
							k = k + 1
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
		DO i = 1,poisson_solver(ps)%real2xpen%ncomm
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
END SUBROUTINE poisson_solver_push_3array
