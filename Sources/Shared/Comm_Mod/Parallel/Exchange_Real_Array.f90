!==============================================================================!
  subroutine Exchange_Real_Array(Comm, length, phi, dest)
!------------------------------------------------------------------------------!
!   Exchanges the values of a real array between the processors.               !
!   (Check out the Sendrecv variant of this function - it is more flexible ... !
!    ... because buffer in each subdomain don't have to have the same length)  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm
  integer,          intent(in)    :: length
  real,             intent(inout) :: phi(length)
  integer,          intent(in)    :: dest         ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer          :: rtag, stag, error  ! receive and send tags, error
  type(Mpi_Status) :: status
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Form send and receive tags
  stag = Communicator % n_processors * Communicator % this_processor + dest
  rtag = Communicator % n_processors * dest + Communicator % this_processor

  call Mpi_Sendrecv_Replace(phi(1),                & ! buffer
                            length,                & ! length
                            comm_type_real,        & ! datatype
                            (dest-1),              & ! dest,
                            stag,                  & ! sendtag,
                            (dest-1),              & ! source,
                            rtag,                  & ! recvtag,
                            MPI_COMM_WORLD,        &
                            status,                &
                            error)

  end subroutine
