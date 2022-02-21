!==============================================================================!
  subroutine Control_Mod_Read_Strings_On(keyword, values, n, verbose)
!------------------------------------------------------------------------------!
!   Used to read variable names in bnd. and initial condition specificaton     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: keyword
  character(SL)     :: values(MSI)   ! spefified value, if found
  integer           :: n             ! number of items
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  logical :: reached_end
!==============================================================================!

  ! Set default values
  values(1:MSI) = ''

  !---------------------------------------------------------!
  !   Read one line from command file to find the keyword   !
  !---------------------------------------------------------!
  call File % Read_Line(control_file_unit, reached_end)
  if(reached_end) goto 1

  ! Found the correct keyword
  if(line % tokens(1) .eq. trim(keyword)) then

    do i=2, line % n_tokens
      read(line % tokens(i), *) values(i-1)
    end do
    n = line % n_tokens - 1
    return 

  ! Keyword not found, try to see if there is similar, maybe it was a typo
  ! (Tokens 2 and on hold variable names, they are too short to be checked)
  else
    call Control_Mod_Similar_Warning( keyword, trim(line % tokens(1)) )
  end if

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(present(verbose)) then
    if(verbose .and. this_proc < 2) then
      print '(a,a,a)', ' # NOTE! Could not find the keyword: ', keyword, '.'
    end if
  end if
  n = 0

  end subroutine
