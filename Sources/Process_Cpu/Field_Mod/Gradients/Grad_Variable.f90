!==============================================================================!
  subroutine Grad_Variable(Flow, var)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to calculate the gradient of a generic
!>  variable in a flow field.  It branches to one of the two methods for
!>  gradient computation: Gauss's theorem-based method and the least squares
!>  method and the choice is based on the grad_method attribute of the generic
!>  variable var.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow  !! parent flow object
  type(Var_Type)    :: var   !! variable whose gradient is calculated
!==============================================================================!

  ! Branch into the method of choice for the variable
  if(var % grad_method .eq. GAUSS_THEOREM) then
    call Flow % Grad_Gauss_Variable(var)
  else
    call Flow % Grad_Least_Variable(var)
  end if

  end subroutine
