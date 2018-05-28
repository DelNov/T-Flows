!==============================================================================!
  subroutine Control_Mod_Max_Simple_Iterations(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_SIMPLE_ITERATIONS', 10,  &
                                  val, verbose)

  end subroutine
