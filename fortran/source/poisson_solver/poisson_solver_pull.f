!------------------------------------------------------------------------------!
!  
!  File:         poisson_solver_pull.f
!  
!  Description:  Pulls the field arrays from the poisson solver object
!  
!------------------------------------------------------------------------------!
SUBROUTINE poisson_solver_pull_1pointer3d(ps,offset,real_lhs,real_rhs)

IMPLICIT NONE

include 'mpif.h'
!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN)                                       :: ps
	INTEGER,DIMENSION(3),INTENT(IN)                            :: offset
	REAL(MK),DIMENSION(:,:,:,:),POINTER,INTENT(INOUT)          :: real_lhs
	REAL(MK),DIMENSION(:,:,:,:),POINTER,INTENT(INOUT),OPTIONAL :: real_rhs

!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER :: nproc, rank, ierr
	INTEGER :: i, j, k
	INTEGER :: p, q, s, sn, sqn, pqs
	INTEGER :: nx, ny
	INTEGER :: iproc, jproc
	INTEGER :: nsend, nrecv
	INTEGER :: nmap
	INTEGER,DIMENSION(3) :: ncell, imin, imax

	LOGICAL :: pack
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
! Check that the pointers are associated
!------------------------------------------------------------------------------!
	IF ( .NOT.ASSOCIATED(real_lhs) ) THEN
		WRITE(*,'(A)')"[poisson_solver_pull]: Error pointer 1 not associated."
		GO TO 9999
	END IF
	IF(PRESENT(real_rhs))THEN
		IF ( .NOT.ASSOCIATED(real_rhs) ) THEN
			WRITE(*,'(A)')"[poisson_solver_pull]: Error pointer 2 not associated."
			GO TO 9999
		END IF
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

	IF ( SIZE(real_lhs,1) .GE. 1 .AND. ALLOCATED(poisson_solver(ps)%lhsX) ) THEN
		lX = .TRUE.
		nmap = nmap + 1
	END IF
	IF ( SIZE(real_lhs,1) .GE. 2 .AND. ALLOCATED(poisson_solver(ps)%lhsY) ) THEN
		lY = .TRUE.
		nmap = nmap + 1
	END IF
	IF ( SIZE(real_lhs,1) .GE. 3 .AND. ALLOCATED(poisson_solver(ps)%lhsZ) ) THEN
		lZ = .TRUE.
		nmap = nmap + 1
	END IF

	IF(PRESENT(real_rhs))THEN
		IF ( SIZE(real_rhs,1) .GE. 1 .AND. ALLOCATED(poisson_solver(ps)%rhsX) ) THEN
			rX = .TRUE.
			nmap = nmap + 1
		END IF
		IF ( SIZE(real_rhs,1) .GE. 2 .AND. ALLOCATED(poisson_solver(ps)%rhsY) ) THEN
			rY = .TRUE.
			nmap = nmap + 1
		END IF
		IF ( SIZE(real_rhs,1) .GE. 3 .AND. ALLOCATED(poisson_solver(ps)%rhsZ) ) THEN
			rZ = .TRUE.
			nmap = nmap + 1
		END IF
	END IF

!------------------------------------------------------------------------------!
! Allocate and pack communication buffers
!------------------------------------------------------------------------------!
	ALLOCATE( buffer_recv(poisson_solver(ps)%xpen2real%ncomm) )

	DO i = 1,poisson_solver(ps)%xpen2real%ncomm

!------------------------------------------------------------------------------!
! Self-communication
!------------------------------------------------------------------------------!
		IF ( poisson_solver(ps)%xpen2real%info(i)%nway .EQ. 0 ) THEN

			imin(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(1)
			imin(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(2)
			imin(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(3)

			imax(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(1)
			imax(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(2)
			imax(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(3)

			nx = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_send(1)
			ny = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_send(2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Fill receive buffer directly
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			nsend = nmap * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(1) &
			             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(2) &
			             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(3)

			ALLOCATE( buffer_recv(i)%val(0:nsend-1) )
			buffer_recv(i)%val = 0.0_MK

			k = 0
			DO  s = imin(3),imax(3)
				sn = s * ny
				DO  q = imin(2),imax(2)
					sqn = ( sn + q ) * nx
					DO  p = imin(1),imax(1)
						pqs = sqn + p

						IF ( rX ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%rhsX(pqs),MK)
							k = k + 1
						END IF
						IF ( rY ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%rhsY(pqs),MK)
							k = k + 1
						END IF
						IF ( rZ ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%rhsZ(pqs),MK)
							k = k + 1
						END IF

						IF ( lX ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%lhsX(pqs),MK)
							k = k + 1
						END IF
						IF ( lY ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%lhsY(pqs),MK)
							k = k + 1
						END IF
						IF ( lZ ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%lhsZ(pqs),MK)
							k = k + 1
						END IF

					END DO
				END DO
			END DO

!------------------------------------------------------------------------------!
! INTEGER ::er-communication
!------------------------------------------------------------------------------!
		ELSE

			iproc = poisson_solver(ps)%xpen2real%info(i)%iproc
			jproc = poisson_solver(ps)%xpen2real%info(i)%jproc

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			pack = .TRUE.
			IF ( rank .EQ. iproc .AND. &
			     poisson_solver(ps)%xpen2real%info(i)%nway .NE. 2 ) THEN

				imin(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(1)
				imin(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(2)
				imin(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(3)

				imax(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(1)
				imax(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(2)
				imax(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(3)

				nsend = nmap * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(1) &
				             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(2) &
				             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(3)

				nx = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_send(1)
				ny = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_send(2)

			ELSE IF ( rank .EQ. jproc .AND. &
			          poisson_solver(ps)%xpen2real%info(i)%nway .GT. 1 ) THEN

				imin(1) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_send(1)
				imin(2) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_send(2)
				imin(3) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_send(3)

				imax(1) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_send(1)
				imax(2) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_send(2)
				imax(3) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_send(3)

				nsend = nmap * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(1) &
				             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(2) &
				             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(3)

				nx = poisson_solver(ps)%xpen2real%info(i)%j2i_ncell_partition_send(1)
				ny = poisson_solver(ps)%xpen2real%info(i)%j2i_ncell_partition_send(2)
			ELSE
				pack = .FALSE.
			END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Pack send buffer
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			IF (pack) THEN
				ALLOCATE( buffer_send%val(0:nsend-1) )
				buffer_send%val = 0.0_MK
				k = 0

				DO  s = imin(3),imax(3)
					sn = s * ny
					DO  q = imin(2),imax(2)
						sqn = ( sn + q ) * nx
						DO  p = imin(1),imax(1)
							pqs = sqn + p

							IF ( rX ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%rhsX(pqs),MK)
								k = k + 1
							END IF
							IF ( rY ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%rhsY(pqs),MK)
								k = k + 1
							END IF
							IF ( rZ ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%rhsZ(pqs),MK)
								k = k + 1
							END IF
							IF ( lX ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%lhsX(pqs),MK)
								k = k + 1
							END IF
							IF ( lY ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%lhsY(pqs),MK)
								k = k + 1
							END IF
							IF ( lZ ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%lhsZ(pqs),MK)
								k = k + 1
							END IF

						END DO
					END DO
				END DO

			END IF

!------------------------------------------------------------------------------!
! Send/Recieve
!------------------------------------------------------------------------------!
			IF ( poisson_solver(ps)%xpen2real%info(i)%nway .EQ. 1 ) THEN ! only i to j
				IF ( rank .EQ. iproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, jproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. jproc ) THEN ! Receiver
					nrecv = nmap * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(1) &
					             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(2) &
					             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, iproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( poisson_solver(ps)%xpen2real%info(i)%nway .EQ. 2 ) THEN ! only j to i
				IF ( rank .EQ. jproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, iproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. iproc ) THEN ! Receiver
					nrecv = nmap * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(1) &
					             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(2) &
					             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( poisson_solver(ps)%xpen2real%info(i)%nway .EQ. 3 ) THEN ! i to j and j to i
				IF ( rank .EQ. iproc ) THEN
					nrecv = nmap * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(1) &
					             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(2) &
					             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_SENDRECV( buffer_send%val, nsend, MPI_DOUBLE, jproc, iproc, &
					              buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, jproc, &
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				ELSE IF ( rank .EQ. jproc ) THEN
					nrecv = nmap * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(1) &
					             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(2) &
					             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(3)
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
! Clear arrays
!------------------------------------------------------------------------------!
	nrecv = poisson_solver(ps)%xpen2real%partition_recv(rank)%ncell(1) &
	      * poisson_solver(ps)%xpen2real%partition_recv(rank)%ncell(2) &
	      * poisson_solver(ps)%xpen2real%partition_recv(rank)%ncell(3)


	DO k = 0,poisson_solver(ps)%xpen2real%partition_recv(rank)%ncell(3)-1
		DO j = 0,poisson_solver(ps)%xpen2real%partition_recv(rank)%ncell(2)-1
			DO i = 0,poisson_solver(ps)%xpen2real%partition_recv(rank)%ncell(1)-1

				IF ( rX ) THEN
					real_rhs(1,i+offset(1),j+offset(2),k+offset(3)) = 0.0_MK
				END IF
				IF ( rY ) THEN
					real_rhs(2,i+offset(1),j+offset(2),k+offset(3)) = 0.0_MK
				END IF
				IF ( rZ ) THEN
					real_rhs(3,i+offset(1),j+offset(2),k+offset(3)) = 0.0_MK
				END IF

				IF ( lX ) THEN
					real_lhs(1,i+offset(1),j+offset(2),k+offset(3)) = 0.0_MK
				END IF
				IF ( lY ) THEN
					real_lhs(2,i+offset(1),j+offset(2),k+offset(3)) = 0.0_MK
				END IF
				IF ( lZ ) THEN
					real_lhs(3,i+offset(1),j+offset(2),k+offset(3)) = 0.0_MK
				END IF

			END DO
		END DO
	END DO

!------------------------------------------------------------------------------!
! Unpack communication buffer
!------------------------------------------------------------------------------!
	DO  i = 1,poisson_solver(ps)%xpen2real%ncomm

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		pack = .TRUE.
		iproc = poisson_solver(ps)%xpen2real%info(i)%iproc
		jproc = poisson_solver(ps)%xpen2real%info(i)%jproc

		IF ( rank .EQ. jproc .AND. &
		     poisson_solver(ps)%xpen2real%info(i)%nway .NE. 2 ) THEN

			imin(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_recv(1)
			imin(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_recv(2)
			imin(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_recv(3)

			imax(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_recv(1)
			imax(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_recv(2)
			imax(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_recv(3)

			nx = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_recv(1)
			ny = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_recv(2)

		ELSE IF ( rank .EQ. iproc .AND. &
		          poisson_solver(ps)%xpen2real%info(i)%nway .GT. 1 ) THEN

			imin(1) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_recv(1)
			imin(2) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_recv(2)
			imin(3) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_recv(3)

			imax(1) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_recv(1)
			imax(2) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_recv(2)
			imax(3) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_recv(3)

			nx = poisson_solver(ps)%xpen2real%info(i)%j2i_ncell_partition_recv(1)
			ny = poisson_solver(ps)%xpen2real%info(i)%j2i_ncell_partition_recv(2)

		ELSE
			pack = .FALSE.
		END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Unpack buffers
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		IF (pack) THEN

			k = 0

			DO  s = imin(3),imax(3)
				sn = s * ny
				DO  q = imin(2),imax(2)
					sqn = ( sn + q ) * nx
					DO  p = imin(1),imax(1)
						pqs = sqn + p

						IF ( rX ) THEN
							real_rhs(1,p+offset(1),q+offset(2),s+offset(3)) = &
							                                          buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( rY ) THEN
							real_rhs(2,p+offset(1),q+offset(2),s+offset(3)) = &
							                                          buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( rZ ) THEN
							real_rhs(3,p+offset(1),q+offset(2),s+offset(3)) = &
							                                          buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( lX ) THEN
							real_lhs(1,p+offset(1),q+offset(2),s+offset(3)) = &
							                                          buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( lY ) THEN
							real_lhs(2,p+offset(1),q+offset(2),s+offset(3)) = &
							                                          buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( lZ ) THEN
							real_lhs(3,p+offset(1),q+offset(2),s+offset(3)) = &
							                                          buffer_recv(i)%val(k)
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
		DO i = 1,poisson_solver(ps)%xpen2real%ncomm
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
END SUBROUTINE poisson_solver_pull_1pointer3d









SUBROUTINE poisson_solver_pull_3array( ps, real_lhsX, real_lhsY, real_lhsZ, &
                                           real_rhsX, real_rhsY, real_rhsZ )

IMPLICIT NONE

include 'mpif.h'

!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	INTEGER,  INTENT(IN) :: ps

	REAL(MK),DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: real_lhsX
	REAL(MK),DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: real_lhsY
	REAL(MK),DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: real_lhsZ

	REAL(MK),DIMENSION(:),ALLOCATABLE,INTENT(INOUT),OPTIONAL :: real_rhsX
	REAL(MK),DIMENSION(:),ALLOCATABLE,INTENT(INOUT),OPTIONAL :: real_rhsY
	REAL(MK),DIMENSION(:),ALLOCATABLE,INTENT(INOUT),OPTIONAL :: real_rhsZ

!------------------------------------------------------------------------------!
! Local variables
!------------------------------------------------------------------------------!
	INTEGER :: nproc, rank, ierr
	INTEGER :: i, j, k
	INTEGER :: p, q, s, sn, sqn, pqs
	INTEGER :: nx, ny
	INTEGER :: iproc, jproc
	INTEGER :: nsend, nrecv
	INTEGER :: nmap
	INTEGER,DIMENSION(3) :: ncell, imin, imax

	LOGICAL :: pack
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
	IF ( ALLOCATED(real_lhsX) .AND. ALLOCATED(poisson_solver(ps)%lhsX) ) THEN
		lX = .TRUE.
		nmap = nmap + 1
	END IF
	IF ( ALLOCATED(real_lhsY) .AND. ALLOCATED(poisson_solver(ps)%lhsY) ) THEN
		lY = .TRUE.
		nmap = nmap + 1
	END IF
	IF ( ALLOCATED(real_lhsZ) .AND. ALLOCATED(poisson_solver(ps)%lhsZ) ) THEN
		lZ = .TRUE.
		nmap = nmap + 1
	END IF

	IF( PRESENT(real_rhsX) )THEN
		IF ( ALLOCATED(real_rhsX) .AND. ALLOCATED(poisson_solver(ps)%rhsX) ) THEN
			rX = .TRUE.
			nmap = nmap + 1
		END IF
	END IF
	IF( PRESENT(real_rhsY) )THEN
		IF ( ALLOCATED(real_rhsY) .AND. ALLOCATED(poisson_solver(ps)%rhsY) ) THEN
			rY = .TRUE.
			nmap = nmap + 1
		END IF
	END IF
	IF( PRESENT(real_rhsZ) )THEN
		IF ( ALLOCATED(real_rhsZ) .AND. ALLOCATED(poisson_solver(ps)%rhsZ) ) THEN
			rZ = .TRUE.
			nmap = nmap + 1
		END IF
	END IF

!------------------------------------------------------------------------------!
! Allocate and pack communication buffers
!------------------------------------------------------------------------------!
	ALLOCATE( buffer_recv(poisson_solver(ps)%xpen2real%ncomm) )

	DO  i = 1,poisson_solver(ps)%xpen2real%ncomm

!------------------------------------------------------------------------------!
! Self-communication
!------------------------------------------------------------------------------!
		IF ( poisson_solver(ps)%xpen2real%info(i)%nway .EQ. 0 ) THEN

			imin(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(1)
			imin(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(2)
			imin(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(3)

			imax(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(1)
			imax(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(2)
			imax(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(3)

			nx = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_send(1)
			ny = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_send(2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Fill receive buffer directly
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			nsend = nmap * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(1) &
			             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(2) &
			             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(3)

			ALLOCATE( buffer_recv(i)%val(0:nsend-1) )
			buffer_recv(i)%val = 0.0_MK

			k = 0
			DO  s = imin(3),imax(3)
				sn = s * ny
				DO  q = imin(2),imax(2)
					sqn = ( sn + q ) * nx
					DO  p = imin(1),imax(1)
						pqs = sqn + p

						IF ( rX ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%rhsX(pqs),MK)
							k = k + 1
						END IF
						IF ( rY ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%rhsY(pqs),MK)
							k = k + 1
						END IF
						IF ( rZ ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%rhsZ(pqs),MK)
							k = k + 1
						END IF

						IF ( lX ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%lhsX(pqs),MK)
							k = k + 1
						END IF
						IF ( lY ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%lhsY(pqs),MK)
							k = k + 1
						END IF
						IF ( lZ ) THEN
							buffer_recv(i)%val(k) = REAL(poisson_solver(ps)%lhsZ(pqs),MK)
							k = k + 1
						END IF

					END DO
				END DO
			END DO

!------------------------------------------------------------------------------!
! INTEGER ::er-communication
!------------------------------------------------------------------------------!
		ELSE

			iproc = poisson_solver(ps)%xpen2real%info(i)%iproc
			jproc = poisson_solver(ps)%xpen2real%info(i)%jproc

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			pack = .TRUE.
			IF ( rank .EQ. iproc .AND. &
			     poisson_solver(ps)%xpen2real%info(i)%nway .NE. 2 ) THEN

				imin(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(1)
				imin(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(2)
				imin(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_send(3)

				imax(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(1)
				imax(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(2)
				imax(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_send(3)

				nsend = nmap * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(1) &
				             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(2) &
				             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(3)

				nx = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_send(1)
				ny = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_send(2)

			ELSE IF ( rank .EQ. jproc .AND. &
			          poisson_solver(ps)%xpen2real%info(i)%nway .GT. 1 ) THEN

				imin(1) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_send(1)
				imin(2) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_send(2)
				imin(3) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_send(3)

				imax(1) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_send(1)
				imax(2) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_send(2)
				imax(3) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_send(3)

				nsend = nmap * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(1) &
				             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(2) &
				             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(3)

				nx = poisson_solver(ps)%xpen2real%info(i)%j2i_ncell_partition_send(1)
				ny = poisson_solver(ps)%xpen2real%info(i)%j2i_ncell_partition_send(2)
			ELSE
				pack = .FALSE.
			END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Pack send buffer
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			IF (pack) THEN
				ALLOCATE( buffer_send%val(0:nsend-1) )
				buffer_send%val = 0.0_MK
				k = 0

				DO  s = imin(3),imax(3)
					sn = s * ny
					DO  q = imin(2),imax(2)
						sqn = ( sn + q ) * nx
						DO  p = imin(1),imax(1)
							pqs = sqn + p

							IF ( rX ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%rhsX(pqs),MK)
								k = k + 1
							END IF
							IF ( rY ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%rhsY(pqs),MK)
								k = k + 1
							END IF
							IF ( rZ ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%rhsZ(pqs),MK)
								k = k + 1
							END IF
							IF ( lX ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%lhsX(pqs),MK)
								k = k + 1
							END IF
							IF ( lY ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%lhsY(pqs),MK)
								k = k + 1
							END IF
							IF ( lZ ) THEN
								buffer_send%val(k) = REAL(poisson_solver(ps)%lhsZ(pqs),MK)
								k = k + 1
							END IF

						END DO
					END DO
				END DO

			END IF

!------------------------------------------------------------------------------!
! Send/Recieve
!------------------------------------------------------------------------------!
			IF ( poisson_solver(ps)%xpen2real%info(i)%nway .EQ. 1 ) THEN ! only i to j
				IF ( rank .EQ. iproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, jproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. jproc ) THEN ! Receiver
					nrecv = nmap * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(1) &
					             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(2) &
					             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, iproc, jproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( poisson_solver(ps)%xpen2real%info(i)%nway .EQ. 2 ) THEN ! only j to i
				IF ( rank .EQ. jproc ) THEN ! Sender
					CALL MPI_SEND( buffer_send%val, nsend, MPI_DOUBLE, iproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				ELSE IF ( rank .EQ. iproc ) THEN ! Receiver
					nrecv = nmap * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(1) &
					             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(2) &
					             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_RECV( buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, iproc, &
					          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				END IF
			ELSE IF ( poisson_solver(ps)%xpen2real%info(i)%nway .EQ. 3 ) THEN ! i to j and j to i
				IF ( rank .EQ. iproc ) THEN
					nrecv = nmap * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(1) &
					             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(2) &
					             * poisson_solver(ps)%xpen2real%info(i)%j2i_ncell(3)
					ALLOCATE( buffer_recv(i)%val( 0:nrecv-1 ) )
					buffer_recv(i)%val = 0.0_MK
					CALL MPI_SENDRECV( buffer_send%val, nsend, MPI_DOUBLE, jproc, iproc, &
					              buffer_recv(i)%val, nrecv, MPI_DOUBLE, jproc, jproc, &
					              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				ELSE IF ( rank .EQ. jproc ) THEN
					nrecv = nmap * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(1) &
					             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(2) &
					             * poisson_solver(ps)%xpen2real%info(i)%i2j_ncell(3)
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
! Clear arrays
!------------------------------------------------------------------------------!
	nrecv = poisson_solver(ps)%xpen2real%partition_recv(rank)%ncell(1) &
	      * poisson_solver(ps)%xpen2real%partition_recv(rank)%ncell(2) &
	      * poisson_solver(ps)%xpen2real%partition_recv(rank)%ncell(3)

	DO i = 0,nrecv-1
		IF ( rX ) THEN
			real_rhsX(i) = 0.0_MK
		END IF
		IF ( rY ) THEN
			real_rhsY(i) = 0.0_MK
		END IF
		IF ( rZ ) THEN
			real_rhsZ(i) = 0.0_MK
		END IF

		IF ( lX ) THEN
			real_lhsX(i) = 0.0_MK
		END IF
		IF ( lY ) THEN
			real_lhsY(i) = 0.0_MK
		END IF
		IF ( lZ ) THEN
			real_lhsZ(i) = 0.0_MK
		END IF
	END DO

!------------------------------------------------------------------------------!
! Unpack communication buffer
!------------------------------------------------------------------------------!
	DO  i = 1,poisson_solver(ps)%xpen2real%ncomm

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Mesh info
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		pack = .TRUE.
		iproc = poisson_solver(ps)%xpen2real%info(i)%iproc
		jproc = poisson_solver(ps)%xpen2real%info(i)%jproc

		IF ( rank .EQ. jproc .AND. &
		     poisson_solver(ps)%xpen2real%info(i)%nway .NE. 2 ) THEN

			imin(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_recv(1)
			imin(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_recv(2)
			imin(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_min_recv(3)

			imax(1) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_recv(1)
			imax(2) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_recv(2)
			imax(3) = poisson_solver(ps)%xpen2real%info(i)%i2j_max_recv(3)

			nx = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_recv(1)
			ny = poisson_solver(ps)%xpen2real%info(i)%i2j_ncell_partition_recv(2)

		ELSE IF ( rank .EQ. iproc .AND. &
		          poisson_solver(ps)%xpen2real%info(i)%nway .GT. 1 ) THEN

			imin(1) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_recv(1)
			imin(2) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_recv(2)
			imin(3) = poisson_solver(ps)%xpen2real%info(i)%j2i_min_recv(3)

			imax(1) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_recv(1)
			imax(2) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_recv(2)
			imax(3) = poisson_solver(ps)%xpen2real%info(i)%j2i_max_recv(3)

			nx = poisson_solver(ps)%xpen2real%info(i)%j2i_ncell_partition_recv(1)
			ny = poisson_solver(ps)%xpen2real%info(i)%j2i_ncell_partition_recv(2)

		ELSE
			pack = .FALSE.
		END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Unpack buffers
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
		IF (pack) THEN

			k = 0

			DO  s = imin(3),imax(3)
				sn = s * ny
				DO  q = imin(2),imax(2)
					sqn = ( sn + q ) * nx
					DO  p = imin(1),imax(1)
						pqs = sqn + p

						IF ( rX ) THEN
							real_rhsX(pqs) = buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( rY ) THEN
							real_rhsY(pqs) = buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( rZ ) THEN
							real_rhsZ(pqs) = buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( lX ) THEN
							real_lhsX(pqs) = buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( lY ) THEN
							real_lhsY(pqs) = buffer_recv(i)%val(k)
							k = k + 1
						END IF
						IF ( lZ ) THEN
							real_lhsZ(pqs) = buffer_recv(i)%val(k)
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
		DO i = 1,poisson_solver(ps)%xpen2real%ncomm
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
END SUBROUTINE poisson_solver_pull_3array
