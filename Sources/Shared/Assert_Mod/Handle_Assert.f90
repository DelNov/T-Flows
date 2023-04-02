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
  character(16) :: numb, proc
!==============================================================================!

  ! Sequential run
  if(Sequential_Run()) then
    if(fail) then
      write(numb, '(i16)') line
      print '(7a)', '  Assert(',         text,                 &
                    ') failed in file ', file,                 &
                    ' at line ',         trim(adjustl(numb)),  &
                    '.'
      call Global % End_Parallel
      stop
    end if

  ! Parallel run
  else
    if(fail) then
      write(numb, '(i16)') line
      write(proc, '(i16)') This_Proc()
      print '(9a)', '  Assert(',         text,                 &
                    ') failed in file ', file,                 &
                    ' at line ',         trim(adjustl(numb)),  &
                    ' in processor ',    trim(adjustl(proc)),  &
                    '.'
      call Global % End_Parallel
      stop
    end if

  end if

  end subroutine
