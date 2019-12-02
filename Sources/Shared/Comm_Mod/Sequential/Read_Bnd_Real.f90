!==============================================================================!
  subroutine Comm_Mod_Read_Bnd_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of reading a "distributed" boundary cell-based array.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type) :: comm
  integer         :: fh         ! file handle
  real            :: array(:)
  integer         :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" boundary cell array 
  do c = 1, comm % nb_t
    read(fh) array(c)
  end do

  disp = disp + comm % nb_t * SIZE_REAL

  end subroutine
