!==============================================================================!
  subroutine Control_Mod_Multiphase_Model(multiphase_model, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: multiphase_model
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('MULTIPHASE_MODEL', 'none', val, verbose)
  call To_Upper_Case(val)

  end subroutine
