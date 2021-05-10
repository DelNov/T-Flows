!==============================================================================!
  subroutine Statistics(Cpu_Timer)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Cpu_Timer_Type), target :: Cpu_Timer
!-----------------------------------[Locals]-----------------------------------!
  integer            :: f
  real               :: total_time, t_temp
  integer            :: hours, minutes, seconds
  character(DL)      :: line, n_temp
  integer, parameter :: I=33                ! indent
  logical            :: swap
!==============================================================================!

  ! Compute average time spent in functions over all processors
  if(n_proc > 1) then
    do f=1, Cpu_Timer % n_funct
      call Comm_Mod_Global_Sum_Real(Cpu_Timer % funct_time(f))
      Cpu_Timer % funct_time(f) = Cpu_Timer % funct_time(f) / n_proc
    end do
  end if

  ! Perform bubble sort
  do
    swap = .false.
    do f=1, Cpu_Timer % n_funct-1
      if(Cpu_Timer % funct_time(f+1) > Cpu_Timer % funct_time(f)) then
        t_temp = Cpu_Timer % funct_time(f)
        n_temp = Cpu_Timer % funct_name(f)
        Cpu_Timer % funct_time(f)   = Cpu_Timer % funct_time(f+1)
        Cpu_Timer % funct_name(f)   = Cpu_Timer % funct_name(f+1)
        Cpu_Timer % funct_time(f+1) = t_temp
        Cpu_Timer % funct_name(f+1) = n_temp
        swap = .true.
      end if
    end do
    if(.not. swap) goto 1
  end do
1 continue

  total_time = 0.0
  do f=1, Cpu_Timer % n_funct
    total_time = total_time + Cpu_Timer % funct_time(f)
  end do

  if(this_proc < 2) then

    line( 1:160) = ' '
    line( 1+I:63+I) =   &
               '#=============================================================#'
    print *, trim(line)
    line( 1+I:63+I) =   &
               '#                    CPU usage statistics                     #'
    print *, trim(line)
    line( 1+I:63+I) =   &
               '#-------------------------------------------------------------#'
    print *, trim(line)
    line( 1:160) = ' '
    line( 1+I:30+I) = '#            Total CPU time:  '

    ! Work out hours, minutes and seconds
    hours   = floor(  total_time  / 3600.0 )
    minutes = floor( (total_time - 3600.0 * hours) / 60.0)
    seconds = floor(  total_time - 3600.0 * hours - 60.0 * minutes )
    write(line(30+I:32+I), '(i3.3)')  hours
    write(line(33+I:33+I),   '(a1)')  ':'
    write(line(34+I:36+I), '(i2.2)')  minutes
    write(line(36+I:36+I),   '(a1)')  ':'
    write(line(37+I:38+I), '(i2.2)')  seconds
    write(line(40+I:50+I),  '(a11)') '[hhh:mm:ss]'
    line(63+I:63+I) = '#'
    print *, trim(line)
    line( 1+I:63+I) =  &
               '#-------------------------------------------+-----------------#'
    print *, trim(line)
    line( 1+I:63+I) =  &
               '#       Description of the activity:        |   Spent time:   #'
    print *, trim(line)
    line( 1+I:63+I) =   &
               '#-------------------------------------------+-----------------#'
    print *, trim(line)

    do f=1, Cpu_Timer % n_funct
      line( 1:160) = ' '
      line( 1+I: 1+I) = '#'
      line(63+I:63+I) = '#'
      line( 3+I: 3+I) = '-'
      line( 5+I: 5+I+len_trim(Cpu_Timer % funct_name(f)))  &
                            = Cpu_Timer % funct_name(f)(1:123)
      line(45+I:45+I) = '|'
      write(line(50+I:55+I), '(f6.2)') Cpu_Timer % funct_time(f)  &
                                       / total_time * 100.0
      line(57+I:57+I) = '%'
      print *, trim(line)
    end do

    line( 1+I:63+I) =  &
               '#-------------------------------------------+-----------------#'
    print *, trim(line)
    print *, ''

  end if

  end subroutine
