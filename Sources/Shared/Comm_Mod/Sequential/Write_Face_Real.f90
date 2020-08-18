!==============================================================================!
  subroutine Comm_Mod_Write_Face_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of writing a "distributed" face-based array.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type) :: comm
  integer         :: fh         ! file handle
  real            :: array(:)
  integer         :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: s
!==============================================================================!

  ! Write "distributed" face data
  do s = 1, comm % nf_t
    write(9) array(s)
  end do

  disp = disp + comm % nf_t * SIZE_REAL

  end subroutine
