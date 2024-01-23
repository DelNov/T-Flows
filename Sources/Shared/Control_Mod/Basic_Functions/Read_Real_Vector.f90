!==============================================================================!
  subroutine Read_Real_Vector(Control, keyword, n, def, val, verbose)
!------------------------------------------------------------------------------!
!   Working horse function to read integer value (argument "val") behind a     !
!   keyword (argument "keyword") in control file.  If not found, a default     !
!   value specified in argument "def" is used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)              :: Control
  character(len=*),    intent(in)  :: keyword
  integer,             intent(in)  :: n        ! size of array (typically small)
  real,                intent(in)  :: def(n)   ! default value
  real,                intent(out) :: val(n)   ! spefified value, if found
  logical,   optional, intent(in)  :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
  integer :: i
!==============================================================================!

  rewind(Control % file_unit)

  ! Set default values
  val = def

  !-----------------------------------------------------!
  !   Browse through command file to find the keyword   !
  !-----------------------------------------------------!
  do
    call File % Read_Line(Control % file_unit, reached_end)
    if(reached_end) goto 1

    ! Found the correct keyword
    if(Line % tokens(1) .eq. trim(keyword)) then
      do i = 1, n
        read(Line % tokens(i+1), *) val(i)
      end do
      return
    end if

  end do

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(present(verbose)) then
    if(verbose .and. First_Proc()) then
      print '(3a,1pe10.3)', ' # NOTE! Could not find the keyword: ',  &
                             trim(keyword), '. Using the default: ', def(1)
    end if
  end if

  end subroutine
