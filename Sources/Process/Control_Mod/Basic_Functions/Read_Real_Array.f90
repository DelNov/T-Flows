!==============================================================================!
  subroutine Control_Mod_Read_Real_Array(keyword, n, def, val, verbose)
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
  integer           :: n        ! size of array (typically small)
  real              :: def(n)   ! default value
  real              :: val(n)   ! spefified value, if found
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
  integer :: i
!==============================================================================!

  rewind(CONTROL_FILE_UNIT)

  ! Set default values
  val = def

  !-----------------------------------------------------!
  !   Browse through command file to find the keyword   !
  !-----------------------------------------------------!
  do
    call Tokenizer_Mod_Read_Line(CONTROL_FILE_UNIT, reached_end)
    if(reached_end) goto 1

    ! Found the correct keyword
    if(line % tokens(1) .eq. trim(keyword)) then
      do i = 1, n
        read(line % tokens(i+1), *) val(i)
      end do
      return

    ! Keyword not found, try to see if there is similar, maybe it was a typo
    else
      call Control_Mod_Similar_Warning( keyword, trim(line % tokens(1)) )
    end if

  end do

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(present(verbose)) then
    if(verbose .and. this_proc < 2) then
      print '(3a,1pe9.3)', ' # NOTE! Could not find the keyword: ',  &
                            trim(keyword), '. Using the default: ', def(1)
    end if
  end if

  end subroutine
