!==============================================================================!
  subroutine Comm_Mod_Write_Bnd_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of writing a "distributed" boundary-cell-based array.   !
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

  ! Write "distributed" boundary cell data 
  do c = 1, comm % nb_tot
    write(9) array(c)
  end do

  disp = disp + comm % nb_tot * SIZE_REAL

  end subroutine
