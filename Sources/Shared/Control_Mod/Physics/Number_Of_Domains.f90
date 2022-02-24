!==============================================================================!
  subroutine Control_Mod_Number_Of_Domains(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('NUMBER_OF_DOMAINS', 1, &
                                  val, verbose)

  end subroutine
