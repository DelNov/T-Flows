!==============================================================================!
  subroutine Comm_Mod_Read_Face_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of reading a "distributed" face-based array.            !
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

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" face data
  do s = 1, comm % nf_t
    read(fh) array(s)
  end do

  disp = disp + comm % nf_t * SIZE_REAL

  end subroutine
