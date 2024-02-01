!==============================================================================!
  character(SL) function Var_Name(str1, str2, bool)
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  character(*) :: str1, str2
  logical      :: bool
!==============================================================================!

  Var_Name = trim(str1)
  if(bool) Var_Name = trim(Var_Name) // " " // trim(str2)

  end function

