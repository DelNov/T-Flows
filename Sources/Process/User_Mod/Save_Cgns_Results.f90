!==============================================================================!
  subroutine User_Mod_Save_Cgns_Results(base, block, solution, field, grid)
!------------------------------------------------------------------------------!
!   This function is called from Save_Cgns and allows user to save his         !
!   own arrays in the Vtu file with results.                                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Var_Mod
  use Cgns_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer         :: base, block, solution, field
  type(Grid_Type) :: grid
!----------------------------------[Locals]------------------------------------!
  integer           :: us, ua
  character(len=80) :: a_name, v_name 
!==============================================================================!

  !-----------------------!
  !   Save user scalars   !
  !-----------------------!
  do us = 1, n_user_scalars
    v_name = 'UserScalarXXXX'
    write(v_name(11:14), '(a4)') user_scalar(us) % name
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              user_scalar(us) % n(1), v_name)
  end do
  
  !----------------------!
  !   Save user arrays   !
  !----------------------!
  do ua = 1, n_user_arrays  
    a_name = 'UserArrayA_XX'
    write(a_name(12:13), '(i2.2)') ua
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              user_array(ua,1), a_name)
  end do

  end subroutine 
