!==============================================================================!
  subroutine Handle_Assert(fail, text, file, line)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical      :: fail
  character(*) :: text
  character(*) :: file
  integer      :: line
!-----------------------------------[Locals]-----------------------------------!
  character(16) :: numb
!==============================================================================!

  if(fail) then
    write(numb, '(i16)') line
    print '(a)', 'Assertion ',       text,                 &
                 ' failed in file ', file,                 &
                 ' at line ',        trim(adjustl(numb)),  &
                 '.'
    call Comm_Mod_End
    stop
  end if

  end subroutine
