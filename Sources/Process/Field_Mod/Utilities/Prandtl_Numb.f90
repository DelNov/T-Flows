!==============================================================================!
  pure real function Prandtl_Numb(Flow, c)
!------------------------------------------------------------------------------!
!   Computes laminar Prandtl number for cell "c"                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(in) :: Flow
  integer,           intent(in) :: c
!==============================================================================!

  Prandtl_Numb = Flow % viscosity(c) * Flow % capacity(c)  &
               / Flow % conductivity(c)

  end function
