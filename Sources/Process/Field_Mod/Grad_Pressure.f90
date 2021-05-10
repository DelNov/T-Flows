!==============================================================================!
  subroutine Field_Mod_Grad_Pressure(flow, p)
!------------------------------------------------------------------------------!
!   Calculates pressure (or pressure correction) gradient field flow           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  type(Var_Type)   :: p
!==============================================================================!

  ! Branch into the method of choice for pressure
  if(p % grad_method .eq. GAUSS_THEOREM) then
    call Field_Mod_Grad_Gauss_Pressure(flow, p)
  else
    call Field_Mod_Grad_Least_Pressure(flow, p)
  end if

  end subroutine
