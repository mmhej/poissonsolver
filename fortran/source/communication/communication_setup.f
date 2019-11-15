!------------------------------------------------------------------------------!
!  
!  File:         setup.f
!  Description:  Sets up the communication class
!  
!------------------------------------------------------------------------------!
SUBROUTINE communication_setup( partition_send, partition_recv, communicator)

USE poisson_solver_partition

IMPLICIT NONE

include 'mpif.h'
!------------------------------------------------------------------------------!
! Arguments
!------------------------------------------------------------------------------!
	TYPE(class_partition), DIMENSION(:),ALLOCATABLE, INTENT(IN):: partition_send
	TYPE(class_partition), DIMENSION(:),ALLOCATABLE, INTENT(IN):: partition_recv
	TYPE(class_communication), INTENT(OUT):: communicator

!------------------------------------------------------------------------------!
! Variables
!------------------------------------------------------------------------------!
	INTEGER :: nproc, rank, ierr
	INTEGER :: i, j, k
	INTEGER :: iproc, jproc
	INTEGER :: icomm, ncomm, sum_comm
	INTEGER, DIMENSION(3) :: imin, imax, jmin, jmax
	INTEGER, DIMENSION(3) :: ncell_send, ncell_recv

	TYPE(class_cominfo):: tmp_info
	TYPE(class_partinfo):: prt_info
	TYPE(class_cominfo),DIMENSION(:),ALLOCATABLE :: tmp_info_all
	TYPE(class_communication) :: tmp_comm

	INTEGER,DIMENSION(:),ALLOCATABLE :: proc_busy
	INTEGER,DIMENSION(:),ALLOCATABLE :: comm_list

!------------------------------------------------------------------------------!
! Get MPI info
!------------------------------------------------------------------------------!
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

!------------------------------------------------------------------------------!
! ALLOCATE
!------------------------------------------------------------------------------!
	ALLOCATE( communicator%partition_send( 0:nproc-1 ) )
	ALLOCATE( communicator%partition_recv( 0:nproc-1 ) )
	ALLOCATE( tmp_info_all(nproc*nproc) )

!------------------------------------------------------------------------------!
! Store partition info
!------------------------------------------------------------------------------!
	DO i = 0,nproc-1
		communicator%partition_send(i)%ncell(1) = partition_send(i)%ncell(1)
		communicator%partition_send(i)%ncell(2) = partition_send(i)%ncell(2)
		communicator%partition_send(i)%ncell(3) = partition_send(i)%ncell(3)

		communicator%partition_recv(i)%ncell(1) = partition_recv(i)%ncell(1)
		communicator%partition_recv(i)%ncell(2) = partition_recv(i)%ncell(2)
		communicator%partition_recv(i)%ncell(3) = partition_recv(i)%ncell(3)
	END DO


!------------------------------------------------------------------------------!
! Find number of communications
!------------------------------------------------------------------------------!
	icomm = 0

	DO iproc = 0, nproc-1
		DO jproc = iproc, nproc-1
 
!------------------------------------------------------------------------------!
! Processors and number of communication ways
!------------------------------------------------------------------------------!
			tmp_info%iproc = iproc
			tmp_info%jproc = jproc
			tmp_info%nway = 0

!------------------------------------------------------------------------------!
! Overlap of subdomains from iproc to jproc
!------------------------------------------------------------------------------!
			imin(1) = partition_send(iproc)%icell(1)
			imin(2) = partition_send(iproc)%icell(2)
			imin(3) = partition_send(iproc)%icell(3)

			imax(1) = partition_send(iproc)%icell(1) + partition_send(iproc)%ncell(1) - 1
			imax(2) = partition_send(iproc)%icell(2) + partition_send(iproc)%ncell(2) - 1
			imax(3) = partition_send(iproc)%icell(3) + partition_send(iproc)%ncell(3) - 1

			jmin(1) = partition_recv(jproc)%icell(1)
			jmin(2) = partition_recv(jproc)%icell(2)
			jmin(3) = partition_recv(jproc)%icell(3)

			jmax(1) = partition_recv(jproc)%icell(1) + partition_recv(jproc)%ncell(1) - 1
			jmax(2) = partition_recv(jproc)%icell(2) + partition_recv(jproc)%ncell(2) - 1
			jmax(3) = partition_recv(jproc)%icell(3) + partition_recv(jproc)%ncell(3) - 1

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Check overlap
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			IF( jmax(1) .GE. imin(1) .AND. jmax(2) .GE. imin(2) .AND. &
			    jmax(3) .GE. imin(3) .AND. imax(1) .GE. jmin(1) .AND. &
			    imax(2) .GE. jmin(2) .AND. imax(3) .GE. jmin(3) ) THEN

				tmp_info%nway = tmp_info%nway + 1

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Overlaping indexes - min index
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				IF ( imin(1) .EQ. jmin(1) ) THEN
					tmp_info%i2j_min_send(1) = 0
					tmp_info%i2j_min_recv(1) = 0
				ELSE IF ( imin(1) .GT. jmin(1) ) THEN
					tmp_info%i2j_min_send(1) = 0
					tmp_info%i2j_min_recv(1) = imin(1) - jmin(1)
				ELSE IF ( imin(1) .LT. jmin(1) ) THEN
					tmp_info%i2j_min_send(1) = jmin(1) - imin(1)
					tmp_info%i2j_min_recv(1) = 0
				END IF

				IF ( imin(2) .EQ. jmin(2) ) THEN
					tmp_info%i2j_min_send(2) = 0
					tmp_info%i2j_min_recv(2) = 0
				ELSE IF ( imin(2) .GT. jmin(2) ) THEN
					tmp_info%i2j_min_send(2) = 0
					tmp_info%i2j_min_recv(2) = imin(2) - jmin(2)
				ELSE IF ( imin(2) .LT. jmin(2) ) THEN
					tmp_info%i2j_min_send(2) = jmin(2) - imin(2)
					tmp_info%i2j_min_recv(2) = 0
				END IF

				IF ( imin(3) .EQ. jmin(3) ) THEN
					tmp_info%i2j_min_send(3) = 0
					tmp_info%i2j_min_recv(3) = 0
				ELSE IF ( imin(3) .GT. jmin(3) ) THEN
					tmp_info%i2j_min_send(3) = 0
					tmp_info%i2j_min_recv(3) = imin(3) - jmin(3)
				ELSE IF ( imin(3) .LT. jmin(3) ) THEN
					tmp_info%i2j_min_send(3) = jmin(3) - imin(3)
					tmp_info%i2j_min_recv(3) = 0
				END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Overlaping indexes - max index
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				IF ( imax(1) .EQ. jmax(1) ) THEN
					tmp_info%i2j_max_send(1) = partition_send(iproc)%ncell(1) - 1
					tmp_info%i2j_max_recv(1) = partition_recv(jproc)%ncell(1) - 1
				ELSE IF ( imax(1) .LT. jmax(1) ) THEN
					tmp_info%i2j_max_send(1) = partition_send(iproc)%ncell(1) - 1
					tmp_info%i2j_max_recv(1) = partition_recv(jproc)%ncell(1) &
					                         - 1 - (jmax(1) - imax(1))
				ELSE IF ( imax(1) .GT. jmax(1) ) THEN
					tmp_info%i2j_max_send(1) = partition_send(iproc)%ncell(1) &
					                         - 1 - (imax(1) - jmax(1))
					tmp_info%i2j_max_recv(1) = partition_recv(jproc)%ncell(1) - 1
				END IF

				IF ( imax(2) .EQ. jmax(2) ) THEN
					tmp_info%i2j_max_send(2) = partition_send(iproc)%ncell(2) - 1
					tmp_info%i2j_max_recv(2) = partition_recv(jproc)%ncell(2) - 1
				ELSE IF ( imax(2) .LT. jmax(2) ) THEN
					tmp_info%i2j_max_send(2) = partition_send(iproc)%ncell(2) - 1
					tmp_info%i2j_max_recv(2) = partition_recv(jproc)%ncell(2) &
					                         - 1 - (jmax(2) - imax(2))
				ELSE IF ( imax(2) .GT. jmax(2) ) THEN
					tmp_info%i2j_max_send(2) = partition_send(iproc)%ncell(2) &
					                         - 1 - (imax(2) - jmax(2))
					tmp_info%i2j_max_recv(2) = partition_recv(jproc)%ncell(2) - 1
				END IF

				IF ( imax(3) .EQ. jmax(3) ) THEN
					tmp_info%i2j_max_send(3) = partition_send(iproc)%ncell(3) - 1
					tmp_info%i2j_max_recv(3) = partition_recv(jproc)%ncell(3) - 1
				ELSE IF ( imax(3) .LT. jmax(3) ) THEN
					tmp_info%i2j_max_send(3) = partition_send(iproc)%ncell(3) - 1
					tmp_info%i2j_max_recv(3) = partition_recv(jproc)%ncell(3) &
					                         - 1 - (jmax(3) - imax(3))
				ELSE IF ( imax(3) .GT. jmax(3) ) THEN
					tmp_info%i2j_max_send(3) = partition_send(iproc)%ncell(3) &
					                         - 1 - (imax(3) - jmax(3))
					tmp_info%i2j_max_recv(3) = partition_recv(jproc)%ncell(3) - 1
				END IF


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store number of cells to be communicated
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				ncell_send(1) = tmp_info%i2j_max_send(1) - tmp_info%i2j_min_send(1) + 1
				ncell_send(2) = tmp_info%i2j_max_send(2) - tmp_info%i2j_min_send(2) + 1
				ncell_send(3) = tmp_info%i2j_max_send(3) - tmp_info%i2j_min_send(3) + 1

				ncell_recv(1) = tmp_info%i2j_max_recv(1) - tmp_info%i2j_min_recv(1) + 1
				ncell_recv(2) = tmp_info%i2j_max_recv(2) - tmp_info%i2j_min_recv(2) + 1
				ncell_recv(3) = tmp_info%i2j_max_recv(3) - tmp_info%i2j_min_recv(3) + 1

				IF ( ncell_send(1) .NE. ncell_recv(1) &
				.OR. ncell_send(2) .NE. ncell_recv(2) &
				.OR. ncell_send(3) .NE. ncell_recv(3) ) THEN
!					std::cout << "ERROR!! number of send and recieve cells does not match" << std::endl
				END IF

				tmp_info%i2j_ncell(1) = ncell_send(1)
				tmp_info%i2j_ncell(2) = ncell_send(2)
				tmp_info%i2j_ncell(3) = ncell_send(3)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store total number of cells of the two partition
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				tmp_info%i2j_ncell_partition_send(1) = partition_send(iproc)%ncell(1)
				tmp_info%i2j_ncell_partition_send(2) = partition_send(iproc)%ncell(2)
				tmp_info%i2j_ncell_partition_send(3) = partition_send(iproc)%ncell(3)

				tmp_info%i2j_ncell_partition_recv(1) = partition_recv(jproc)%ncell(1)
				tmp_info%i2j_ncell_partition_recv(2) = partition_recv(jproc)%ncell(2)
				tmp_info%i2j_ncell_partition_recv(3) = partition_recv(jproc)%ncell(3)

			END IF

!------------------------------------------------------------------------------!
! Overlap of subdomains from jproc to iproc
!------------------------------------------------------------------------------!
			jmin(1) = partition_send(jproc)%icell(1)
			jmin(2) = partition_send(jproc)%icell(2)
			jmin(3) = partition_send(jproc)%icell(3)

			jmax(1) = partition_send(jproc)%icell(1) + partition_send(jproc)%ncell(1) - 1
			jmax(2) = partition_send(jproc)%icell(2) + partition_send(jproc)%ncell(2) - 1
			jmax(3) = partition_send(jproc)%icell(3) + partition_send(jproc)%ncell(3) - 1

			imin(1) = partition_recv(iproc)%icell(1)
			imin(2) = partition_recv(iproc)%icell(2)
			imin(3) = partition_recv(iproc)%icell(3)

			imax(1) = partition_recv(iproc)%icell(1) + partition_recv(iproc)%ncell(1) - 1
			imax(2) = partition_recv(iproc)%icell(2) + partition_recv(iproc)%ncell(2) - 1
			imax(3) = partition_recv(iproc)%icell(3) + partition_recv(iproc)%ncell(3) - 1



!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Check overlap
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
			IF( imax(1) .GE. jmin(1) .AND. imax(2) .GE. jmin(2) .AND. &
			    imax(3) .GE. jmin(3) .AND. jmax(1) .GE. imin(1) .AND. &
			    jmax(2) .GE. imin(2) .AND. jmax(3) .GE. imin(3) .AND. &
			    iproc .NE. jproc ) THEN

				tmp_info%nway = tmp_info%nway + 2

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Overlaping indexes - min index
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				IF ( jmin(1) .EQ. imin(1) ) THEN
					tmp_info%j2i_min_send(1) = 0
					tmp_info%j2i_min_recv(1) = 0
				ELSE IF ( jmin(1) .GT. imin(1) ) THEN
					tmp_info%j2i_min_send(1) = 0
					tmp_info%j2i_min_recv(1) = jmin(1) - imin(1)
				ELSE IF ( jmin(1) .LT. imin(1) ) THEN
					tmp_info%j2i_min_send(1) = imin(1) - jmin(1)
					tmp_info%j2i_min_recv(1) = 0
				END IF

				IF ( jmin(2) .EQ. imin(2) ) THEN
					tmp_info%j2i_min_send(2) = 0
					tmp_info%j2i_min_recv(2) = 0
				ELSE IF ( jmin(2) .GT. imin(2) ) THEN
					tmp_info%j2i_min_send(2) = 0
					tmp_info%j2i_min_recv(2) = jmin(2) - imin(2)
				ELSE IF ( jmin(2) .LT. imin(2) ) THEN
					tmp_info%j2i_min_send(2) = imin(2) - jmin(2)
					tmp_info%j2i_min_recv(2) = 0
				END IF

				IF ( jmin(3) .EQ. imin(3) ) THEN
					tmp_info%j2i_min_send(3) = 0
					tmp_info%j2i_min_recv(3) = 0
				ELSE IF ( jmin(3) .GT. imin(3) ) THEN
					tmp_info%j2i_min_send(3) = 0
					tmp_info%j2i_min_recv(3) = jmin(3) - imin(3)
				ELSE IF ( jmin(3) .LT. imin(3) ) THEN
					tmp_info%j2i_min_send(3) = imin(3) - jmin(3)
					tmp_info%j2i_min_recv(3) = 0
				END IF


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Overlaping indexes - max index
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				IF ( jmax(1) .EQ. imax(1) ) THEN
					tmp_info%j2i_max_send(1) = partition_send(jproc)%ncell(1) - 1
					tmp_info%j2i_max_recv(1) = partition_recv(iproc)%ncell(1) - 1
				ELSE IF ( jmax(1) .LT. imax(1) ) THEN
					tmp_info%j2i_max_send(1) = partition_send(jproc)%ncell(1) - 1
					tmp_info%j2i_max_recv(1) = partition_recv(iproc)%ncell(1) &
					                         - 1 - (imax(1) - jmax(1))
				ELSE IF ( jmax(1) .GT. imax(1) ) THEN
					tmp_info%j2i_max_send(1) = partition_send(jproc)%ncell(1) &
					                         - 1 - (jmax(1) - imax(1))
					tmp_info%j2i_max_recv(1) = partition_recv(iproc)%ncell(1) - 1
				END IF

				IF ( jmax(2) .EQ. imax(2) ) THEN
					tmp_info%j2i_max_send(2) = partition_send(jproc)%ncell(2) - 1
					tmp_info%j2i_max_recv(2) = partition_recv(iproc)%ncell(2) - 1
				ELSE IF ( jmax(2) .LT. imax(2) ) THEN
					tmp_info%j2i_max_send(2) = partition_send(jproc)%ncell(2) - 1
					tmp_info%j2i_max_recv(2) = partition_recv(iproc)%ncell(2) &
					                         - 1 - (imax(2) - jmax(2))
				ELSE IF ( jmax(2) .GT. imax(2) ) THEN
					tmp_info%j2i_max_send(2) = partition_send(jproc)%ncell(2) &
					                         - 1 - (jmax(2) - imax(2))
					tmp_info%j2i_max_recv(2) = partition_recv(iproc)%ncell(2) - 1
				END IF

				IF ( jmax(3) .EQ. imax(3) ) THEN
					tmp_info%j2i_max_send(3) = partition_send(jproc)%ncell(3) - 1
					tmp_info%j2i_max_recv(3) = partition_recv(iproc)%ncell(3) - 1
				ELSE IF ( jmax(3) .LT. imax(3) ) THEN
					tmp_info%j2i_max_send(3) = partition_send(jproc)%ncell(3) - 1
					tmp_info%j2i_max_recv(3) = partition_recv(iproc)%ncell(3) &
					                         - 1 - (imax(3) - jmax(3))
				ELSE IF ( jmax(3) .GT. imax(3) ) THEN
					tmp_info%j2i_max_send(3) = partition_send(jproc)%ncell(3) &
					                         - 1 - (jmax(3) - imax(3))
					tmp_info%j2i_max_recv(3) = partition_recv(iproc)%ncell(3) - 1
				END IF



!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store number of cells to be communicated
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				ncell_send(1) = tmp_info%j2i_max_send(1) - tmp_info%j2i_min_send(1) + 1
				ncell_send(2) = tmp_info%j2i_max_send(2) - tmp_info%j2i_min_send(2) + 1
				ncell_send(3) = tmp_info%j2i_max_send(3) - tmp_info%j2i_min_send(3) + 1

				ncell_recv(1) = tmp_info%j2i_max_recv(1) - tmp_info%j2i_min_recv(1) + 1
				ncell_recv(2) = tmp_info%j2i_max_recv(2) - tmp_info%j2i_min_recv(2) + 1
				ncell_recv(3) = tmp_info%j2i_max_recv(3) - tmp_info%j2i_min_recv(3) + 1

				IF ( ncell_send(1) .NE. ncell_recv(1) &
				.OR. ncell_send(2) .NE. ncell_recv(2) &
				.OR. ncell_send(3) .NE. ncell_recv(3) ) THEN
!					std::cout << "ERROR!! number of send and recieve cells does not match" << std::endl
				END IF

				tmp_info%j2i_ncell(1) = ncell_send(1)
				tmp_info%j2i_ncell(2) = ncell_send(2)
				tmp_info%j2i_ncell(3) = ncell_send(3)


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Store total number of cells of the two partition
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
				tmp_info%j2i_ncell_partition_send(1) = partition_send(jproc)%ncell(1)
				tmp_info%j2i_ncell_partition_send(2) = partition_send(jproc)%ncell(2)
				tmp_info%j2i_ncell_partition_send(3) = partition_send(jproc)%ncell(3)

				tmp_info%j2i_ncell_partition_recv(1) = partition_recv(iproc)%ncell(1)
				tmp_info%j2i_ncell_partition_recv(2) = partition_recv(iproc)%ncell(2)
				tmp_info%j2i_ncell_partition_recv(3) = partition_recv(iproc)%ncell(3)

!				IF( rank .EQ. 0 .AND. ncell_send(1) .NE. ncell_recv(1) ) THEN
!				IF( rank .EQ. 0 ) THEN
!					WRITE(*,*) "ERROR!! "
!					WRITE(*,*) "proc: ",jproc,iproc
!					WRITE(*,*) "imin: ",jmin(1),imin(1), tmp_info%j2i_min_send(1), tmp_info%j2i_min_recv(1)
!					WRITE(*,*) "imax: ",jmax(1),imax(1), tmp_info%j2i_max_send(1), tmp_info%j2i_max_recv(1)
!					WRITE(*,*) "ncells: ", ncell_send(1),ncell_recv(1)
!				END IF

			END IF


!------------------------------------------------------------------------------!
! Store communication
!------------------------------------------------------------------------------!
			IF (tmp_info%nway .GT. 0) THEN
				IF ( iproc .EQ. jproc ) THEN
					tmp_info%nway = 0
				END IF
				icomm = icomm + 1
				tmp_info_all(icomm) = tmp_info
			END IF

		END DO !jproc
	END DO !iproc

! Store total number of communications
	ncomm = icomm

!------------------------------------------------------------------------------!
! Find communication sequence (a greedy colouring algorithm)
!------------------------------------------------------------------------------!
	ALLOCATE( tmp_comm%info(ncomm) )
	ALLOCATE( proc_busy( 0:nproc-1 ) )
	proc_busy = 0
	ALLOCATE( comm_list(ncomm) )
	comm_list = 1

	tmp_comm%ncomm = 0
	icomm = 0

	sum_comm = ncomm
	DO WHILE ( sum_comm .GT. 0 )

! Free processors
		DO i = 0,nproc-1
			proc_busy(i) = 0
		END DO

! Loop through all communications
		sum_comm = 0
		DO i = 1,ncomm

			iproc = tmp_info_all(i)%iproc
			jproc = tmp_info_all(i)%jproc

! Check if communication is still listed and if processors are available
			IF( comm_list(i) .EQ. 1 .AND. &
			    proc_busy(iproc) .EQ. 0 .AND. proc_busy(jproc) .EQ. 0 ) THEN
				comm_list(i) = 0
				proc_busy(iproc) = 1
				proc_busy(jproc) = 1


!	IF (rank .EQ. 0) THEN
!		WRITE(*,*) "Procs: ", iproc, jproc, tmp_info_all(i)%nway
!	END IF


! Append communication to the processor specific communication list
				IF ( rank .EQ. iproc .OR. rank .EQ. jproc ) THEN
					icomm = icomm + 1
					tmp_comm%info(icomm) = tmp_info_all(i)
					tmp_comm%ncomm = tmp_comm%ncomm + 1
				END IF

			END IF

! Sum remaining communications
			sum_comm = sum_comm + comm_list(i)

		END DO

	END DO

!------------------------------------------------------------------------------!
! Store communication sequence
!------------------------------------------------------------------------------!
	communicator%ncomm = tmp_comm%ncomm
	ALLOCATE( communicator%info(tmp_comm%ncomm) )
	DO i = 1,tmp_comm%ncomm
		communicator%info(i) = tmp_comm%info(i)
	END DO

!------------------------------------------------------------------------------!
! De-allocate
!------------------------------------------------------------------------------!
	DEALLOCATE( tmp_comm%info )
	DEALLOCATE( tmp_info_all )
	DEALLOCATE( proc_busy )
	DEALLOCATE( comm_list )

!------------------------------------------------------------------------------!
! Output communication
!------------------------------------------------------------------------------!
!DO j = 0,nproc
!	IF (rank .EQ. j) THEN
!		WRITE(*,*) "rank: ", rank
!		WRITE(*,*) "total ncoms: ", ncomm
!		WRITE(*,*) "ncom rounds: ", communicator%ncomm
!i = 1
!		DO i = 1,communicator%ncomm
!			WRITE(*,*) "imin: ",communicator%info(i)%i2j_min_send(1), communicator%info(i)%i2j_min_recv(1) &
!			                   ,communicator%info(i)%i2j_min_send(2), communicator%info(i)%i2j_min_recv(2) &
!			                   ,communicator%info(i)%i2j_min_send(3), communicator%info(i)%i2j_min_recv(3)
!			WRITE(*,*) "imax: ",communicator%info(i)%i2j_max_send(1), communicator%info(i)%i2j_max_recv(1) &
!			                   ,communicator%info(i)%i2j_max_send(2), communicator%info(i)%i2j_max_recv(2) &
!			                   ,communicator%info(i)%i2j_max_send(3), communicator%info(i)%i2j_max_recv(3)
!			WRITE(*,*) "--"
!		END DO
!	END IF
!	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!END DO
!					WRITE(*,*) "ncells: ", ncell_send(1),ncell_recv(1)

!------------------------------------------------------------------------------!
! Return
!------------------------------------------------------------------------------!
	9999 CONTINUE
RETURN
END SUBROUTINE communication_setup

