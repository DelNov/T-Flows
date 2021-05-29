!==============================================================================!
  subroutine Recv_Log_Array(Comm, len_r, phi_r, dest)
!------------------------------------------------------------------------------!
!   Receives a logical array from processor dest.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: len_r         ! receive length
  logical          :: phi_r(len_r)  ! receive buffer
  integer          :: dest          ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer :: rtag, error
  integer :: status(MPI_STATUS_SIZE)
!==============================================================================!

  ! Form Receive tags
  rtag = (n_proc) * dest + this_proc  ! tag for receiving

  call Mpi_Recv(phi_r(1),        & ! receive buffer
                len_r,           & ! receive length
                MPI_LOGICAL8,    & ! datatype
                (dest-1),        & ! source,
                rtag,            & ! recvtag,
                MPI_COMM_WORLD,  &
                status,          &
                error)

  end subroutine
