!==============================================================================!
  real function Roughness_Coefficient(grid, z_o_function, c)
!------------------------------------------------------------------------------!
!   Set lower limit to roughness coefficient based on wall distance.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,       only: Grid_Type
  use Rans_Mod,       only: z_o
  use Const_Mod,      only: TINY
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: c
  real            :: z_o_function
!==============================================================================!

  Roughness_Coefficient = z_o

  if(z_o_function > tiny) then
    Roughness_Coefficient = z_o_function
  end if

  end function
