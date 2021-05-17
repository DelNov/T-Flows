!==============================================================================!
  real function Prandtl_Number(Flow, c)
!------------------------------------------------------------------------------!
!   Computes laminar Prandtl number for cell "c"                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  integer                   :: c
!==============================================================================!

  Prandtl_Number = Flow % viscosity(c) * Flow % capacity(c)  &
                 / Flow % conductivity(c)

  end function
