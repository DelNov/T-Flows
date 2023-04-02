!==============================================================================!
  subroutine Control_Mod_Number_Of_Scalars(val, verbose)
!------------------------------------------------------------------------------!
!   Reading stuff related to passive scalars                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_SCALARS', 0, val, verbose)

  end subroutine
