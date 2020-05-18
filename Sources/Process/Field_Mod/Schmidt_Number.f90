!==============================================================================!
  real function Field_Mod_Schmidt_Number(flow, c)
!------------------------------------------------------------------------------!
!   Computes laminar Schmidt number for cell "c"                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  integer                  :: c
!==============================================================================!

  Field_Mod_Schmidt_Number = flow % viscosity(c)  &
                           / (flow % diffusivity * flow % density(c))

  end function
