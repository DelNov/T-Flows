!==============================================================================!
  subroutine Comm_Mod_Send_Int_Array(len_s, phi_s, dest)
!------------------------------------------------------------------------------!
!   Sends an integer array to processor dest.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: len_s         ! send length
  integer :: phi_s(len_s)  ! send buffer
  integer :: dest          ! destination processor
!-----------------------------------[Locals]-----------------------------------!
  integer :: stag, error
  integer :: status(MPI_STATUS_SIZE)
!==============================================================================!

  ! Form send tag
  stag = (n_proc) * this_proc + dest  ! tag for sending

  call Mpi_Send(phi_s(1),        & ! send buffer
                len_s,           & ! send length
                MPI_INTEGER8,    & ! datatype
                (dest-1),        & ! dest,
                stag,            & ! sendtag,
                MPI_COMM_WORLD,  &
                status,          &
                error)

  end subroutine
