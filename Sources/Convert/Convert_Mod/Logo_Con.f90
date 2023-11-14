!==============================================================================!
  subroutine Logo_Con(Convert)
!------------------------------------------------------------------------------!
!>  Prints Convert's logo on the terminal.  Along with it, also in which
!>  precision (single or double) were the floating point numbers compiled
!>  and the flavour of sorting routines used (quick or heap sort).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

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
  print *,'#-----------------------------------' // &
          '------------------------------------'

  end subroutine
