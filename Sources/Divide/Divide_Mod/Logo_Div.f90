!==============================================================================!
  subroutine Logo_Div(Divide)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!---------------------------------[Arguments]----------------------------------!
  class(Divide_Type) :: Divide
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
  print *,'#                         ___  _      _    __'
  print *,'#                        / _ \(_)  __(_)__/ /__'
  print *,'#                       / // / / |/ / / _  / -_)'
  print *,'#                      /____/_/|___/_/\_,_/\__/'
  print *,'#'
  if(RP .eq. DP) then
  print *,'#                        Double precision mode'
  else
  print *,'#                        Single precision mode'
  end if
  print *,'#-----------------------------------' // &
          '------------------------------------'

  end subroutine
