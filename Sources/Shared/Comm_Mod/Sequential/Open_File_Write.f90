!==============================================================================!
  subroutine Comm_Mod_Open_File_Write(fh, file_name)
!------------------------------------------------------------------------------!
!   Open file for writing for sequential runs.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: fh             ! file handle
  character :: file_name*(*)  ! file name
!==============================================================================!

  fh = 9

  open(fh, FILE=file_name, access='stream', form='unformatted')

  end subroutine
