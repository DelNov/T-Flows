!==============================================================================!
  subroutine Exchange_Log_Array(Comm, length, phi, dest)
!------------------------------------------------------------------------------!
!>  This subroutine manages the transfer of logical arrays between processors
!>  in parallel computing, synchronizing data across various subdomains.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm         !! communication object
  integer,          intent(in)    :: length       !! length of the array
  logical,          intent(inout) :: phi(length)  !! array to be exchanged
  integer,          intent(in)    :: dest         !! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer          :: rtag, stag, error  ! receive and send tags, error
  type(Mpi_Status) :: status
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Form send and receive tags
  stag = Global % n_processors * Global % this_processor + dest
  rtag = Global % n_processors * dest + Global % this_processor

  call Mpi_Sendrecv_Replace(phi(1),          & ! buffer
                            length,          & ! length
                            comm_type_log,   & ! datatype
                            (dest-1),        & ! dest,
                            stag,            & ! sendtag,
                            (dest-1),        & ! source,
                            rtag,            & ! recvtag,
                            MPI_COMM_WORLD,  &
                            status,          &
                            error)

  end subroutine
