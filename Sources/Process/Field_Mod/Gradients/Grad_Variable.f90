!==============================================================================!
  subroutine Grad_Variable(Flow, var)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable from field Flow                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
  type(Var_Type)    :: var
!==============================================================================!

  ! Branch into the method of choice for the variable
  if(var % grad_method .eq. GAUSS_THEOREM) then
    call Flow % Grad_Gauss_Variable(var)
  else
    call Flow % Grad_Least_Variable(var)
  end if

  end subroutine
