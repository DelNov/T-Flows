!==============================================================================!
  subroutine Grad_Pressure_Correction(Flow, pp)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure correction.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
  type(Var_Type)    :: pp
!==============================================================================!

  ! Branch into the method of choice for pressure
  if(pp % grad_method .eq. GAUSS_THEOREM) then
    call Flow % Grad_Gauss_Pressure(pp)
  else
    call Flow % Grad_Least_Pressure_Correction(pp)
  end if

  end subroutine
