!==============================================================================!
  subroutine Switch_To_Domain(Control, dom)
!------------------------------------------------------------------------------!
!   Switch the control file to the specified domain.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  integer, intent(in) :: dom
!==============================================================================!

  Control % file_unit = Control % dom_file_unit(dom)

  end subroutine
