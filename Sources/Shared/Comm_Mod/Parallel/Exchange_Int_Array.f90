!==============================================================================!
  subroutine Comm_Mod_Exchange_Int_Array(phi, length, dest)
!------------------------------------------------------------------------------!
!   Exchanges the values of an integer array between the processors.           !
!   (Check out the Sendrecv variant of this function - it is more flexible ... !
!    ... because buffer in each subdomain don't have to have the same length)  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: phi(length)
  integer :: length
  integer :: dest         ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer :: rtag, stag, error
  integer :: status(MPI_STATUS_SIZE)
!==============================================================================!

  ! Form send and receive tags
  stag = (n_proc) * this_proc + dest  ! tag for sending
  rtag = (n_proc) * dest + this_proc  ! tag for receiving

  call Mpi_Sendrecv_Replace(phi(1),          & ! buffer
                            length,          & ! length
                            MPI_INTEGER8,    & ! datatype
                            (dest-1),        & ! dest,
                            stag,            & ! sendtag,
                            (dest-1),        & ! source,
                            rtag,            & ! recvtag,
                            MPI_COMM_WORLD,  &
                            status,          &
                            error)

  end subroutine
