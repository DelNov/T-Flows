!==============================================================================!
  subroutine User_Mod_Allocate(grid)
!------------------------------------------------------------------------------!
!   This function allocates memory for user-defined variables.                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: us, ua
!==============================================================================!

  !-----------------------------!
  !   User scalars and arrays   !
  !-----------------------------!
  call Control_Mod_Number_Of_User_Arrays(grid % n_user_arrays,  &
                                         verbose = .true.)

  !-------------------------------------!
  !   Allocate memory for user arrays   !
  !-------------------------------------!
  allocate(grid % user_array(grid % n_user_arrays,  &
                            -grid % n_bnd_cells:grid % n_cells))
  do ua = 1, grid % n_user_arrays
    grid % user_array(ua,:) = 0.
  end do

  end subroutine
