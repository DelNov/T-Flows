!==============================================================================!
  subroutine Read_Strings_On(Control, keyword, values, n, verbose)
!------------------------------------------------------------------------------!
!   Used to read variable names in bnd. and initial condition specificaton     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  character(len=*)    :: keyword
  character(SL)       :: values(MSI)   ! spefified value, if found
  integer             :: n             ! number of items
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  logical :: reached_end
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Control)
!==============================================================================!

  ! Set default values
  values(1:MSI) = ''

  !---------------------------------------------------------!
  !   Read one line from command file to find the keyword   !
  !---------------------------------------------------------!
  call File % Read_Line(control_file_unit, reached_end)
  if(reached_end) goto 1

  ! Found the correct keyword
  if(Line % tokens(1) .eq. trim(keyword)) then

    do i=2, Line % n_tokens
      read(Line % tokens(i), *) values(i-1)
    end do
    n = Line % n_tokens - 1
    return
  end if

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(present(verbose)) then
    if(verbose .and. First_Proc()) then
      print '(a,a,a)', ' # NOTE! Could not find the keyword: ', keyword, '.'
    end if
  end if
  n = 0

  end subroutine
