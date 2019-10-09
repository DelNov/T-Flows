!==============================================================================!
  subroutine Cgns_Mod_Close_File
!------------------------------------------------------------------------------!
!   Closes file_id file [sequential vesion]                                    !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Close a CGNS file
  call Cg_Close_F(file_id,  & !(in )
                  error)      !(out)

  if (error .ne. 0) then
    print *, '# Failed to close the file: ', trim(file_name)
    call Cg_Error_Exit_F()
  endif

  end subroutine
