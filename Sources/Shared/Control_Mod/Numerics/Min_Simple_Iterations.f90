!==============================================================================!
  subroutine Control_Mod_Min_Simple_Iterations(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MIN_SIMPLE_ITERATIONS', 3,  &
                                  val, verbose)

  end subroutine
