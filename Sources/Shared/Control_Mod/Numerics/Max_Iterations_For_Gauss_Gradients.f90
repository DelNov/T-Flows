!==============================================================================!
  subroutine Control_Mod_Max_Iterations_For_Gauss_Gradients(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_ITERATIONS_FOR_GAUSS_GRADIENTS',  &
                                  val, val, verbose)

  end subroutine
