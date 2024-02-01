!==============================================================================!
  character(SL) function Var_Name(Results, v_name, v_unit, add_unit)
!------------------------------------------------------------------------------!
!>  This functions sets name of a variable to be saved in .vtu file format.
!>  Depending on parameter add_unit, it will either return just the variable
!>  name, or the variable name and its unit.
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Results_Type) :: Results   !! parent class
  character(*)        :: v_name    !! variable name
  character(*)        :: v_unit    !! variable unit
  logical             :: add_unit  !! switch to use unit
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  Var_Name = trim(v_name)
  if(add_unit) Var_Name = trim(Var_Name) // " " // trim(v_unit)

  end function

