!==============================================================================!
  subroutine Logo_Con(Convert)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert
!==============================================================================!

  print *,'#===================================' // &
          '===================================='
  print *,'#                                                                   '
  print *,'#    ______________________ ____    ________  __      __  _________ '
  print *,'#    \__    ___/\_   _____/|    |   \_____  \/  \    /  \/   _____/ '
  print *,'#      |    |    |    __)  |    |    /   |   \   \/\/   /\_____  \  '
  print *,'#      |    |    |     \   |    |___/    |    \        / /        \ '
  print *,'#      |____|    \___  /   |_______ \_______  /\__/\  / /_______  / '
  print *,'#                    \/            \/       \/      \/          \/  '
  print *,'#                     _____                      __                 '
  print *,'#                    / ___/__  ___ _  _____ ____/ /_                '
  print *,'#                   / /__/ _ \/ _ \ |/ / -_) __/ __/                '
  print *,'#                   \___/\___/_//_/___/\__/_/  \__/                 '
  print *,'#                                                                   '
  if(RP .eq. DP) then
  print *,'#                         Double precision mode'
  else
  print *,'#                         Single precision mode'
  end if
  print *,'#-----------------------------------' // &
          '------------------------------------'

  end subroutine