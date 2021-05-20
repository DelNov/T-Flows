!==============================================================================!
  subroutine Open_File_Write(Comm, fh, file_name)
!------------------------------------------------------------------------------!
!   Open file for writing for sequential runs.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh             ! file handle
  character        :: file_name*(*)  ! file name
!==============================================================================!

  fh = 9

  open(fh, FILE=file_name, access='stream', form='unformatted')

  end subroutine
