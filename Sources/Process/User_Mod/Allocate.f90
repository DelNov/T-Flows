!==============================================================================!
  subroutine User_Mod_Allocate(grid)
!------------------------------------------------------------------------------!
!   This function allocates memory for user-defined variables.                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)    :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: us, ua
!==============================================================================!

  !--------------------------------------!
  !   Allocate memory for user scalars   !
  !--------------------------------------!
  allocate(user_scalar(n_user_scalars))

  !-------------------------------------!
  !   Browse through all user scalars   !
  !-------------------------------------!
  do us = 1, n_user_scalars

    ! Set variable name 
    write(c_name(3:4),'(i2.2)') us

    ! Allocate memory for user scalar
    call Var_Mod_Allocate_Solution(c_name, user_scalar(us), grid)

  end do

  !-------------------------------------!
  !   Allocate memory for user arrays   !
  !-------------------------------------!
  allocate(user_array(n_user_arrays, -grid % n_bnd_cells:grid % n_cells))
  do ua = 1, n_user_arrays 
    user_array(ua,:) = 0.
  end do

  end subroutine
