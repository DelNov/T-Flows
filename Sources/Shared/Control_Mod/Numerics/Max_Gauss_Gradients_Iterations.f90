!==============================================================================!
  subroutine Control_Mod_Max_Gauss_Gradients_Iterations(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_GAUSS_GRADIENTS_ITERATIONS', 12,  &
                                val, verbose)

  end subroutine
