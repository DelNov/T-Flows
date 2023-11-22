!==============================================================================!
  subroutine Purge_Sub(File)
!------------------------------------------------------------------------------!
!   Called from Divide only to purge the "Sub" directory structure             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)           :: File
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: rel_path, sys_comm
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  rel_path = 'Sub/'
  sys_comm = 'rm -fR ' // trim(rel_path)
  call system(trim(sys_comm))

  end subroutine
