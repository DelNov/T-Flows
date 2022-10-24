!==============================================================================!
  subroutine Statistics(Profiler, indent)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Profiler_Type), target :: Profiler
  integer, intent(in)          :: indent     ! 33 for Main_Pro, 1 for Main_con
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i_fun
  real(DP)      :: total_time, t_temp, percent_time
  integer       :: hours, minutes, seconds
  character(DL) :: line, n_temp
  logical       :: swap
!==============================================================================!

  ! Compute average time spent in functions over all processors
  if(n_proc > 1) then
    do i_fun=1, Profiler % n_functions
      call Comm_Mod_Global_Sum_Real(Profiler % funct_time(i_fun))
      Profiler % funct_time(i_fun) = Profiler % funct_time(i_fun) / n_proc
    end do
  end if

  ! Perform bubble sort
  do
    swap = .false.
    do i_fun=1, Profiler % n_functions-1
      if(Profiler % funct_time(i_fun+1) > Profiler % funct_time(i_fun)) then
        t_temp = Profiler % funct_time(i_fun)
        n_temp = Profiler % funct_name(i_fun)
        Profiler % funct_time(i_fun)   = Profiler % funct_time(i_fun+1)
        Profiler % funct_name(i_fun)   = Profiler % funct_name(i_fun+1)
        Profiler % funct_time(i_fun+1) = t_temp
        Profiler % funct_name(i_fun+1) = n_temp
        swap = .true.
      end if
    end do
    if(.not. swap) goto 1
  end do
1 continue

  total_time = 0.0
  do i_fun = 1, Profiler % n_functions
    total_time = total_time + Profiler % funct_time(i_fun)
  end do

  if(this_proc < 2) then

    line( 1:160) = ' '
    line( 1+indent:63+indent) =   &
               '#=============================================================#'
    print '(a)', trim(line)
    line( 1+indent:63+indent) =   &
               '#                    CPU usage statistics                     #'
    print '(a)', trim(line)
    line( 1+indent:63+indent) =   &
               '#-------------------------------------------------------------#'
    print '(a)', trim(line)
    line( 1:160) = ' '
    line( 1+indent:30+indent) = '#            Total CPU time:  '

    ! Work out hours, minutes and seconds
    hours   = floor(  total_time  / 3600.0 )
    minutes = floor( (total_time - 3600.0 * hours) / 60.0)
    seconds = floor(  total_time - 3600.0 * hours - 60.0 * minutes )
    write(line(30+indent:32+indent), '(i3.3)')  hours
    write(line(33+indent:33+indent),   '(a1)')  ':'
    write(line(34+indent:36+indent), '(i2.2)')  minutes
    write(line(36+indent:36+indent),   '(a1)')  ':'
    write(line(37+indent:38+indent), '(i2.2)')  seconds
    write(line(40+indent:50+indent),  '(a11)') '[hhh:mm:ss]'
    line(63+indent:63+indent) = '#'
    print '(a)', trim(line)
    line( 1+indent:63+indent) =  &
               '#-------------------------------------------+-----------------#'
    print '(a)', trim(line)
    line( 1+indent:63+indent) =  &
               '#       Description of the activity:        |   Spent time:   #'
    print '(a)', trim(line)
    line( 1+indent:63+indent) =   &
               '#-------------------------------------------+-----------------#'
    print '(a)', trim(line)

    do i_fun = 1, Profiler % n_functions
      line( 1:160) = ' '
      line( 1+indent: 1+indent) = '#'
      line(63+indent:63+indent) = '#'
      line( 3+indent: 3+indent) = '-'
      line( 5+indent: 5+indent+len_trim(Profiler % funct_name(i_fun)))  &
                                = Profiler % funct_name(i_fun)(1:123)
      line(45+indent:45+indent) = '|'
      percent_time = Profiler % funct_time(i_fun) / total_time * 100.0
      write(line(50+indent:55+indent), '(f6.2)') percent_time
      line(57+indent:57+indent) = '%'
      if(percent_time > 0.01) then
        print '(a)', trim(line)
      end if
    end do

    line( 1+indent:63+indent) =  &
               '#-------------------------------------------+-----------------#'
    print '(a)', trim(line)
    print *, ''

  end if

  end subroutine
