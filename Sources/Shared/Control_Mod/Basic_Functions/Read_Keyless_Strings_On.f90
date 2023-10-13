!==============================================================================!
  subroutine Read_Keyless_Strings_On(Control, values, n, verbose)
!------------------------------------------------------------------------------!
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)              :: Control
  character(SL),       intent(out) :: values(MAX_STRING_ITEMS)
  integer,             intent(out) :: n        ! number of items
  logical,   optional, intent(in)  :: verbose
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  logical :: reached_end
!==============================================================================!

  ! Set default values
  values(1:MAX_STRING_ITEMS) = ''

  !-------------------------------------!
  !   Read one line from command file   !
  !-------------------------------------!
  call File % Read_Line(Control % file_unit, reached_end)
  if(reached_end) goto 1

  ! Just fetch the strings
  do i=1, Line % n_tokens
    read(Line % tokens(i), *) values(i)
  end do
  n = Line % n_tokens
  return

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(present(verbose)) then
    if(verbose .and. First_Proc()) then
      print '(a)', ' # NOTE! Could not find the keyless strings'
    end if
  end if
  n = 0

  end subroutine
