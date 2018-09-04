!==============================================================================!
  subroutine Control_Mod_Max_Iterations_For_User_Scalars_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_ITERATIONS_FOR_USER_SCALARS_SOLVER',  &
                                  val, val, verbose)

  end subroutine
