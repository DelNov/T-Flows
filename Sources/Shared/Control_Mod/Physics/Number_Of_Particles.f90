!==============================================================================!
  subroutine Control_Mod_Number_Of_Particles(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('NUMBER_OF_PARTICLES', 0, &
                                  val, verbose)

  end subroutine
