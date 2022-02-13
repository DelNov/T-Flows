!==============================================================================!
  subroutine Exchange_Int_Array(Comm, length, phi, dest)
!------------------------------------------------------------------------------!
!   Exchanges the values of an integer array between the processors.           !
!   (Check out the Sendrecv variant of this function - it is more flexible ... !
!    ... because buffer in each subdomain don't have to have the same length)  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: length
  integer          :: phi(length)
  integer          :: dest         ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer          :: rtag, stag, error
  type(Mpi_Status) :: status
!==============================================================================!

  ! Form send and receive tags
  stag = (n_proc) * this_proc + dest  ! tag for sending
  rtag = (n_proc) * dest + this_proc  ! tag for receiving

  call Mpi_Sendrecv_Replace(phi(1),          & ! buffer
                            length,          & ! length
                            MPI_INTEGER,     & ! datatype
                            (dest-1),        & ! dest,
                            stag,            & ! sendtag,
                            (dest-1),        & ! source,
                            rtag,            & ! recvtag,
                            MPI_COMM_WORLD,  &
                            status,          &
                            error)

  end subroutine
