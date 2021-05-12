!==============================================================================!
  subroutine Control_Mod_Max_Least_Squares_Gradients_Iterations(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_LEAST_SQUARES_GRADIENTS_ITERATIONS', &
                                  val, val, verbose)

  end subroutine
