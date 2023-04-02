!==============================================================================!
  subroutine Switch_To_Domain(Control, dom)
!------------------------------------------------------------------------------!
!   Switch the control file to the specified domain.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  integer, intent(in) :: dom
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Control)
!==============================================================================!

  control_file_unit = dom_control_file_unit(dom)

  end subroutine
