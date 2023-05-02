!==============================================================================!
  subroutine Open_File_Read(Comm, fh, file_name)
!------------------------------------------------------------------------------!
!   Open file for reading for sequential runs.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh             ! file handle
  character        :: file_name*(*)  ! file name
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  fh = 9

  open(fh, FILE=file_name, access='stream', form='unformatted')

  end subroutine
