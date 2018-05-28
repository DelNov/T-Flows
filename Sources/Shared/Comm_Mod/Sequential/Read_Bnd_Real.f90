!==============================================================================!
  subroutine Comm_Mod_Read_Bnd_Real(fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of reading a "distributed" boundary cell-based array.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh         ! file handle
  real    :: array(:)
  integer :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" boundary cell array 
  do c = 1, nb_t
    read(fh) array(c)
  end do

  disp = disp + nb_t * SIZE_REAL

  end subroutine
