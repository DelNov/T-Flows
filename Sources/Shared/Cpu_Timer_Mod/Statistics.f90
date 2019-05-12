!==============================================================================!
  subroutine Cpu_Timer_Mod_Statistics
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer            :: f
  real               :: total_time, t_temp
  character(len=160) :: line, n_temp
  integer, parameter :: I=33                ! indent
  logical            :: swap
!==============================================================================!

  ! Perform bubble sort
  do
    swap = .false.
    do f=1, n_funct-1
      if(funct_time(f+1) > funct_time(f)) then
        t_temp = funct_time(f)
        n_temp = funct_name(f)
        funct_time(f)   = funct_time(f+1)
        funct_name(f)   = funct_name(f+1)
        funct_time(f+1) = t_temp
        funct_name(f+1) = n_temp
        swap = .true.
      end if
    end do
    if(.not. swap) goto 1
  end do
1 continue

  total_time = 0.0
  do f=1, n_funct
    total_time = total_time + funct_time(f)
  end do

  line( 1:160) = " "
  line( 1+I:63+I) =   &
               "#=============================================================#"
  print *, trim(line)
  line( 1+I:63+I) =   &
               "#                    CPU usage statistics                     #"
  print *, trim(line)
  line( 1+I:63+I) =   &
               "#-------------------------------------------------------------#"
  print *, trim(line)
  line( 1:160) = " "
  line( 1+I:32+I) = "#                Total CPU time: "
  write(line(30+I:39+I), "(f9.3)") total_time
  line(40+I:42+I) = "[s]"
  line(63+I:63+I) = "#"
  print *, trim(line)
  line( 1+I:63+I) =  &
               "#-------------------------------------------------------------#"
  print *, trim(line)
  line( 1+I:63+I) =  &
               "# Function:                               Time:               #"
  print *, trim(line)
  line( 1+I:63+I) =   &
               "#-------------------------------------------------------------#"
  print *, trim(line)

  do f=1, n_funct
    line( 1:160) = " "
    line( 1+I: 1+I) = "#"
    line(63+I:63+I) = "#"
    line( 3+I: 3+I+len_trim(funct_name(f))) = funct_name(f)
    line(39+I:39+I) = ":"
    write(line(41+I:51+I), "(f9.3)") funct_time(f)
    line(51+I:53+I) = "[s]"
    write(line(55+I:59+I), "(f5.1)") funct_time(f) / total_time * 100.0
    line(61+I:61+I) = "%"
    print *, trim(line)
  end do

  line( 1+I:63+I) =  &
               "#-------------------------------------------------------------#"
  print *, trim(line)
  print *, ""

  end subroutine
