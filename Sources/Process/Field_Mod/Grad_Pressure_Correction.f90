!==============================================================================!
  subroutine Field_Mod_Grad_Pressure_Correction(flow, pp)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure correction.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  type(Var_Type)   :: pp
!==============================================================================!

  ! Branch into the method of choice for pressure
  if(pp % grad_method .eq. GAUSS_THEOREM) then
    call Field_Mod_Grad_Gauss_Pressure(flow, pp)
  else
    call Field_Mod_Grad_Least_Pressure_Correction(flow, pp)
  end if

  end subroutine
