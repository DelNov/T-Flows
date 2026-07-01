!==============================================================================!
  subroutine Logo_Con(Convert)
!------------------------------------------------------------------------------!
!>  Prints Convert's logo on the terminal.  Along with it, it reports
!>  the floating-point precision used at compile time, the sorting
!>  routine selected at compile time, and the date and time at which
!>  Convert was launched.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert
!-----------------------------------[Locals]-----------------------------------!
  character(8)  :: date
  character(10) :: time
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  call date_and_time(date=date, time=time)

  print *,'#===================================' // &
          '===================================='
  print *,'#'
  print *,'#    ______________________ ____    ________  __      __  _________'
  print *,'#    \__    ___/\_   _____/|    |   \_____  \/  \    /  \/   _____/'
  print *,'#      |    |    |    __)  |    |    /   |   \   \/\/   /\_____  \'
  print *,'#      |    |    |     \   |    |___/    |    \        / /        \'
  print *,'#      |____|    \___  /   |_______ \_______  /\__/\  / /_______  /'
  print *,'#                    \/            \/       \/      \/          \/'
  print *,'#                     _____                      __'
  print *,'#                    / ___/__  ___ _  _____ ____/ /_'
  print *,'#                   / /__/ _ \/ _ \ |/ / -_) __/ __/'
  print *,'#                   \___/\___/_//_/___/\__/_/  \__/'
  print *,'#'
  if(RP .eq. DP) then
  print *,'#                         Double precision mode'
  else
  print *,'#                         Single precision mode'
  end if
# if T_FLOWS_QUICKSORT == 1
  print *,'#                  Compiled with recursive quicksort'
# else
  print *,'#                 Compiled with nonrecursive heapsort'
# endif
  print *,'#' // YELLOW //'                   Launched on ' //  &
          date(1:4) // '-' // date(5:6) // '-' // date(7:8) //  &
          ' '                                               //  &
          time(1:2) // ':' // time(3:4) // ':' // time(5:6) //  &
          RESET
  print *,'#-----------------------------------' // &
          '------------------------------------'

  end subroutine
