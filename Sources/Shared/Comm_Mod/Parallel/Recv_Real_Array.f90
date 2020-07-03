!==============================================================================!
  subroutine Comm_Mod_Recv_Real_Array(phi_r, len_r, dest)
!------------------------------------------------------------------------------!
!   Receives a real array from processor dest.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real    :: phi_r(len_r)  ! receive buffer
  integer :: len_r         ! receive length
  integer :: dest          ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer :: rtag, error
  integer :: status(MPI_STATUS_SIZE)
!==============================================================================!

  ! Form Receive tags
  rtag = (n_proc) * dest + this_proc  ! tag for receiving

  call Mpi_Recv(phi_r(1),              & ! receive buffer
                len_r,                 & ! receive length
                MPI_DOUBLE_PRECISION,  & ! datatype
                (dest-1),              & ! source,
                rtag,                  & ! recvtag,
                MPI_COMM_WORLD,        &
                status,                &
                error)

  end subroutine
