!==============================================================================!
  subroutine User_Mod_Allocate(grid)
!------------------------------------------------------------------------------!
!   This function allocates memory for user-defined variables.                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,   only: Grid_Type
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)    :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: us, ua
!==============================================================================!

  !-----------------------------!
  !   User scalars and arrays   !
  !-----------------------------!
  call Control_Mod_Number_Of_User_Arrays (n_user_arrays,  verbose = .true.)

  !-------------------------------------!
  !   Allocate memory for user arrays   !
  !-------------------------------------!
  allocate(user_array(n_user_arrays, -grid % n_bnd_cells:grid % n_cells))
  do ua = 1, n_user_arrays
    user_array(ua,:) = 0.
  end do

  end subroutine
