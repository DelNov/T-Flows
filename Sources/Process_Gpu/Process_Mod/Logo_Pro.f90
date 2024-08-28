!==============================================================================!
  subroutine Logo_Pro(Process)
!------------------------------------------------------------------------------!
!>  Prints Process's logo on the terminal.  Along with it, also prints in which
!>  precision (single or double) were the floating point numbers compiled
!>  and if the Process was compiled with or without PETSc solvers.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type) :: Process
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  if(First_Proc()) then
  print *,'#===================================' // &
          '===================================='
  print *,'#'

  !----------------------------------------------------------------------!
  !   Nvidia compiler can't seem to handle backslash in the output ...   !
  !----------------------------------------------------------------------!
# ifdef __NVCOMPILER
  print *,'#    ______________________ ____    ________  __      __  _________'
  print *,'#    \\__    ___/\\_   _____/|    |   \\_____  \\/  \\    /  \\/   _____/'
  print *,'#      |    |    |    __)  |    |    /   |   \\   \\/\\/   /\\_____  \\'
  print *,'#      |    |    |     \\   |    |___/    |    \\        / /        \\'
  print *,'#      |____|    \\___  /   |_______ \\_______  /\\__/\\  / /_______  /'
  print *,'#                    \\/            \\/       \\/      \\/          \\/'
  print *,'#                      ___'
  print *,'#                     / _ \\_______  _______ ___ ___'
  print *,'#                    / ___/ __/ _ \\/ __/ -_|_-<(_-<'
  print *,'#                   /_/  /_/  \\___/\\__/\\__/___/___/'
  !-------------------------------------!
  !   ... whereas other compilers can   !
  !-------------------------------------!
# else
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
# endif
  print *,'#'
# if T_FLOWS_GPU == 0
  print *,'#                          Compiled for CPUs'
# else
  print *,'#                          Compiled for GPUs'
# endif
  if(RP .eq. DP) then
  print *,'#                        Double precision mode'
  else
  print *,'#                        Single precision mode'
  end if
  print *,'#                    Only native solvers available'
  print *,'#-----------------------------------' // &
          '------------------------------------'
  end if

  end subroutine
