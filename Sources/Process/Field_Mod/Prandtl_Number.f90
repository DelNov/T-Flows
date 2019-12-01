!==============================================================================!
  real function Field_Mod_Prandtl_Number(flow, c)
!------------------------------------------------------------------------------!
!   Computes laminar Prandtl number for cell "c"                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  integer                  :: c
!==============================================================================!

  Field_Mod_Prandtl_Number = flow % viscosity(c) * flow % capacity(c)  &
                           / flow % conductivity(c)

  end function
