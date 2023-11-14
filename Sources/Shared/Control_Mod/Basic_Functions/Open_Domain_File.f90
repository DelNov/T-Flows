!==============================================================================!
  subroutine Open_Domain_File(Control, dom, file_name)
!------------------------------------------------------------------------------!
!>  Opens control file for a specified domain.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)             :: Control    !! parent class
  integer,             intent(in) :: dom        !! domain rank (number)
  character(len=*),    intent(in) :: file_name  !! file name (it should be
                                                !! control.1 ... control.n)
!==============================================================================!

  call File % Open_For_Reading_Ascii(file_name,                     &
                                     Control % dom_file_unit(dom))

  end subroutine
