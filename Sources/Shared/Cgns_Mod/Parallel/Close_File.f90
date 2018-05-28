!==============================================================================!
  subroutine Cgns_Mod_Close_File
!------------------------------------------------------------------------------!
!   Closes file_id file [parallel vesion]                                      !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Close a CGNS file
  call Cgp_Close_F(file_id,  & !(in )
                   error)      !(out)

  if (error .ne. 0) then
    print *, "# Failed to close the file: ", trim(file_name)
    call Cgp_Error_Exit_F()
  endif

  end subroutine