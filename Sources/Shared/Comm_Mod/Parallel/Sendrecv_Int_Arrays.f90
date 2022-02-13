!==============================================================================!
  subroutine Sendrecv_Int_Arrays(Comm, len_s, phi_s,  &
                                       len_r, phi_r, dest)
!------------------------------------------------------------------------------!
!   Sends and receives values of two integer arrays between the processors.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: len_s         ! send length
  integer          :: phi_s(len_s)  ! send buffer
  integer          :: len_r         ! receive length
  integer          :: phi_r(len_r)  ! receive buffer
  integer          :: dest          ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer          :: rtag, stag, error
  type(Mpi_Status) :: status
!==============================================================================!

  ! Form send and receive tags
  stag = (n_proc) * this_proc + dest  ! tag for sending
  rtag = (n_proc) * dest + this_proc  ! tag for receiving

  call Mpi_Sendrecv(phi_s(1),              & ! send buffer
                    len_s,                 & ! send length
                    MPI_INTEGER,           & ! datatype
                    (dest-1),              & ! dest,
                    stag,                  & ! sendtag,
                    phi_r(1),              & ! receive buffer
                    len_r,                 & ! receive length
                    MPI_INTEGER,           & ! datatype
                    (dest-1),              & ! source,
                    rtag,                  & ! recvtag,
                    MPI_COMM_WORLD,        &
                    status,                &
                    error)

  end subroutine
