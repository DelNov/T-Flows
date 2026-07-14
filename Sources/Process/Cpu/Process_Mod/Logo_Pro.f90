!==============================================================================!
  subroutine Logo_Pro(Process)
!------------------------------------------------------------------------------!
!>  Prints Process's logo on the terminal, including compile-time
!>  information about floating-point precision and PETSc solver support,
!>  as well as the date and time at which Process was launched.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type) :: Process
!-----------------------------------[Locals]-----------------------------------!
  character(8)  :: date
  character(10) :: time
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  if(First_Proc()) then

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
  print *,'#                      ___'
  print *,'#                     / _ \_______  _______ ___ ___'
  print *,'#                    / ___/ __/ _ \/ __/ -_|_-<(_-<'
  print *,'#                   /_/  /_/  \___/\__/\__/___/___/'
  print *,'#'
  if(RP .eq. DP) then
  print *,'#                        Double precision mode'
  else
  print *,'#                        Single precision mode'
  end if
  if(PETSC_ACTIVE) then
  print *,'#                     Compiled with PETSc solvers'
  else
  print *,'#                    Only native solvers available'
  end if
  print *,'#' // YELLOW //'                   Launched on ' //  &
          date(1:4) // '-' // date(5:6) // '-' // date(7:8) //  &
          ' '                                               //  &
          time(1:2) // ':' // time(3:4) // ':' // time(5:6) //  &
          RESET
  print *,'#-----------------------------------' // &
          '------------------------------------'
  end if

  end subroutine
