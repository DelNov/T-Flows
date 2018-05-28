!==============================================================================!
  subroutine Control_Mod_Read_Strings_On(keyword, values, n, verbose)
!------------------------------------------------------------------------------!
!   Working horse function to read integer value (argument "val") behind a     !
!   keyword (argument "keyword") in control file.  If not found, a default     !
!   vaue specified in argument "def" is used.
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Comm_Mod, only: this_proc
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: keyword
  character(len=80) :: values(128)   ! spefified value, if found
  integer           :: n             ! number of items
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  logical :: reached_end
!==============================================================================!

  ! Set default values
  values = ''

  ! Read one line only to see if you get expected output
  call Tokenizer_Mod_Read_Line(CONTROL_FILE_UNIT, reached_end)
  if(reached_end) goto 1
  if(line % tokens(1) .eq. trim(keyword)) then

    do i=2, line % n_tokens
      read(line % tokens(i), *) values(i-1)
    end do
    n = line % n_tokens - 1
    return 
  end if

1 if(present(verbose)) then
    if(verbose .and. this_proc < 2) then
      print '(a,a,a)', ' # Could not find the keyword: ', keyword, '.'
    end if
  end if
  n = 0

  end subroutine