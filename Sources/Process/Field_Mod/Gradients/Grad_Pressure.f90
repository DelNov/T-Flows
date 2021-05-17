!==============================================================================!
  subroutine Grad_Pressure(Flow, p)
!------------------------------------------------------------------------------!
!   Calculates pressure (or pressure correction) gradient field flow           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
  type(Var_Type)    :: p
!==============================================================================!

  ! Branch into the method of choice for pressure
  if(p % grad_method .eq. GAUSS_THEOREM) then
    call Flow % Grad_Gauss_Pressure(p)
  else
    call Flow % Grad_Least_Pressure(p)
  end if

  end subroutine
