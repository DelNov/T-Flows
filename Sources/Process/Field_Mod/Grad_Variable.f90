!==============================================================================!
  subroutine Field_Mod_Grad_Variable(flow, var)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable from field flow                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  type(Var_Type)   :: var
!==============================================================================!

  ! Branch into the method of choice for the variable
  if(var % grad_method .eq. GAUSS_THEOREM) then
    call Field_Mod_Grad_Gauss_Variable(flow, var)
  else
    call Field_Mod_Grad_Least_Variable(flow, var)
  end if

  end subroutine
