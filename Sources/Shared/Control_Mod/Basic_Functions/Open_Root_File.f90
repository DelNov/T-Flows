!==============================================================================!
  subroutine Open_Root_File(Control, file_name)
!------------------------------------------------------------------------------!
!>  Opens the root control file.  If it doesn't exist, throws an error.
!>  T-Flows can't work withut the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control    !! parent class
  character(len=*)    :: file_name  !! control file name (in essence: control)
!-----------------------------------[Locals]-----------------------------------!
  logical :: file_exists
!==============================================================================!

  ! First check if the file exists
  inquire(file  = trim(file_name),  &
          exist = file_exists)

  ! File doesn't exist
  if(.not. file_exists) then
    call Message % Error(60, "The control file is not "  //   &
                             " present in the working "  //   &
                             " directory!  \n \n      "  //   &
                             " This error is critical."  //   &
                             " Exiting!",                     &
                             one_proc = .true.)
  end if

  call File % Open_For_Reading_Ascii(file_name, Control % root_file_unit)

  ! Make root default to begin with
  Control % file_unit = Control % root_file_unit

  ! Set default values for domain 1
  Control % dom_file_unit(1) = Control % root_file_unit

  end subroutine
