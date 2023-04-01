!==============================================================================!
  pure real function Schmidt_Numb(Flow, c)
!------------------------------------------------------------------------------!
!   Computes laminar Schmidt number for cell "c"                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(in) :: Flow
  integer,           intent(in) :: c
!==============================================================================!

  Schmidt_Numb =  Flow % viscosity(c)  &
               / (Flow % diffusivity * Flow % density(c))

  end function
