!==============================================================================!
  subroutine Grad_Pressure(Flow, p)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to calculate the gradient of pressure or
!>  pressure correction in a flow field.  It branches to one of the two methods
!>  for gradient computation: Gauss's theorem-based method and the least
!>  squares method and the choice is based on the grad_method attribute of
!>  the p (pressure or pressure correction) variable.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow  !! parent flow object
  type(Var_Type)    :: p     !! pressure or pressure correction
!==============================================================================!

  ! Branch into the method of choice for pressure.  (Pressure gradients
  ! balanced with body forces are formulated for least squares only.)
  if(p % grad_method .eq. GAUSS_THEOREM .and.  &
     .not. (Flow % balanced_p_grad .and. p % name .eq. 'P')) then
    call Flow % Grad_Gauss_Pressure(p)
  else
    call Flow % Grad_Least_Pressure(p)
  end if

  end subroutine
