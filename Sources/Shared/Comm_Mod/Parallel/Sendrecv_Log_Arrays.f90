!==============================================================================!
  subroutine Sendrecv_Log_Arrays(Comm, len_s, phi_s,  &
                                       len_r, phi_r, dest)
!------------------------------------------------------------------------------!
!   Sends and receives values of two logical arrays between the processors.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)  :: Comm
  integer,          intent(in)  :: len_s         ! send length
  logical,          intent(in)  :: phi_s(len_s)  ! send buffer
  integer,          intent(in)  :: len_r         ! receive length
  logical,          intent(out) :: phi_r(len_r)  ! receive buffer
  integer,          intent(in)  :: dest          ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer          :: rtag, stag, error  ! tags to send and receive, error
  type(Mpi_Status) :: status
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Form send and receive tags
  stag = Global % n_processors * Global % this_processor + dest
  rtag = Global % n_processors * dest + Global % this_processor

  call Mpi_Sendrecv(phi_s(1),              & ! send buffer
                    len_s,                 & ! send length
                    comm_type_log,         & ! datatype
                    (dest-1),              & ! dest,
                    stag,                  & ! sendtag,
                    phi_r(1),              & ! receive buffer
                    len_r,                 & ! receive length
                    comm_type_log,         & ! datatype
                    (dest-1),              & ! source,
                    rtag,                  & ! recvtag,
                    MPI_COMM_WORLD,        &
                    status,                &
                    error)

  end subroutine
