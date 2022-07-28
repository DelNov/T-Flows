!==============================================================================!
  subroutine Control_Mod_Max_Simple_Iterations(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_SIMPLE_ITERATIONS', 12,  &
                                  val, verbose)

  end subroutine
