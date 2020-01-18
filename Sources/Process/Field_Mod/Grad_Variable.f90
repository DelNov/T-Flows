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

  call Field_Mod_Grad_Component(flow, var % n, 1, var % x)  ! dp/dx
  call Field_Mod_Grad_Component(flow, var % n, 2, var % y)  ! dp/dy
  call Field_Mod_Grad_Component(flow, var % n, 3, var % z)  ! dp/dz

  end subroutine
