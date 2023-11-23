!==============================================================================!
  subroutine Read_Real_Vector_On(Control, keyword, values, n, verbose)
!------------------------------------------------------------------------------!
!>  Working horse function to read real values (argument "values") behind a
!>  keyword (argument "keyword") in control file, starting from the current
!>  position. If not found, it sets them all to zero, which is hard-coded.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)              :: Control      !! parent class
  character(len=*),    intent(in)  :: keyword      !! keyword it searches
  real,                intent(out) :: values(128)  !! spefified value, if found
  integer,             intent(out) :: n            !! number of items
  logical,   optional, intent(in)  :: verbose      !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
  integer :: i
!==============================================================================!

  ! Set default values
  values = 0.0

  !---------------------------------------------------------!
  !   Read one line from command file to find the keyword   !
  !---------------------------------------------------------!
  call File % Read_Line(Control % file_unit, reached_end)
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
      print '(2a)', ' # NOTE! Could not find the keyword: ',  &
                      trim(keyword)
    end if
  end if
  n = 0

  end subroutine
