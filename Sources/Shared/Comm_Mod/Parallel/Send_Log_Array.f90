!==============================================================================!
  subroutine Comm_Mod_Send_Log_Array(phi_s, len_s, dest)
!------------------------------------------------------------------------------!
!   Sends a logical array to processor dest.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical :: phi_s(len_s)  ! send buffer
  integer :: len_s         ! send length
  integer :: dest          ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer :: stag, error
  integer :: status(MPI_STATUS_SIZE)
!==============================================================================!

  ! Form send tag
  stag = (n_proc) * this_proc + dest  ! tag for sending

  call Mpi_Send(phi_s(1),        & ! send buffer
                len_s,           & ! send length
                MPI_LOGICAL,     & ! datatype
                (dest-1),        & ! dest,
                stag,            & ! sendtag,
                MPI_COMM_WORLD,  &
                status,          &
                error)

  end subroutine
