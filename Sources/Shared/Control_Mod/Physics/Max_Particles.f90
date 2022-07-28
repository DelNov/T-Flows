!==============================================================================!
  subroutine Control_Mod_Max_Particles(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_PARTICLES', 0,  &
                                  val, verbose)

  end subroutine
