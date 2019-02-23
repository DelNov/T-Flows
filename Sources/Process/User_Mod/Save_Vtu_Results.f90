!==============================================================================!
  subroutine User_Mod_Save_Vtu_Results(grid)
!------------------------------------------------------------------------------!
!   This function is called from Save_Vtu and allows user to save his          !
!   own arrays in the Vtu file with results.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
!----------------------------------[Locals]------------------------------------!
  integer          :: us, ua
  character(len=4) :: a_name
!-----------------------------[Local parameters]-------------------------------!
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!==============================================================================!

  !----------------------!
  !   Save user arrays   !
  !----------------------!
  do ua = 1, n_user_arrays

    a_name = 'A_00'
    write(a_name(3:4), '(i2.2)') ua

    call Save_Vtu_Scalar(grid, IN_4, IN_5, a_name, user_array(ua,1))
  end do

  end subroutine  ! fourth level comments
