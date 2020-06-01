!==============================================================================!
  subroutine Comm_Mod_Sendrecv_Log_Arrays(phi_s, len_s,  &
                                          phi_r, len_r, dest)
!------------------------------------------------------------------------------!
!   Sends and receives values of two logical arrays between the processors.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical :: phi_s(len_s)  ! send buffer
  integer :: len_s         ! send length
  logical :: phi_r(len_r)  ! receive buffer
  integer :: len_r         ! receive length
  integer :: dest          ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer :: rtag, stag, error
  integer :: status(MPI_STATUS_SIZE)
!==============================================================================!

  ! Form send and receive tags
  stag = (n_proc) * this_proc + dest  ! tag for sending
  rtag = (n_proc) * dest + this_proc  ! tag for receiving

  call Mpi_Sendrecv(phi_s(1),              & ! send buffer
                    len_s,                 & ! send length
                    MPI_LOGICAL8,          & ! datatype
                    (dest-1),              & ! dest,
                    stag,                  & ! sendtag,
                    phi_r(1),              & ! receive buffer
                    len_r,                 & ! receive length
                    MPI_LOGICAL8,          & ! datatype
                    (dest-1),              & ! source,
                    rtag,                  & ! recvtag,
                    MPI_COMM_WORLD,        &
                    status,                &
                    error)

  end subroutine
