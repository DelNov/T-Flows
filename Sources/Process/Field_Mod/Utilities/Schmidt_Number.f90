!==============================================================================!
  real function Schmidt_Number(Flow, c)
!------------------------------------------------------------------------------!
!   Computes laminar Schmidt number for cell "c"                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  integer                   :: c
!==============================================================================!

  Schmidt_Number =  Flow % viscosity(c)  &
                 / (Flow % diffusivity * Flow % density(c))

  end function
