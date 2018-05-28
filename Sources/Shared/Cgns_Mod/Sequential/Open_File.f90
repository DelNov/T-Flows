!==============================================================================!
  subroutine Cgns_Mod_Open_File(file, mode)
!------------------------------------------------------------------------------!
!   Opens name_in file and set its file_id [sequential vesion]                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*) :: file
  integer :: mode ! CG_MODE_READ, CG_MODE_WRITE or CG_MODE_MODIFY
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  file_name = file

  print *, "# Opening file:", trim(file_name)

  ! Open a CGNS file
  call Cg_Open_F(file_name,  & !(in )
                 mode,       & !(in )
                 file_id,    & !(out)
                 error)        !(out)

  if (error .ne. 0) then
    print *, "# Failed to open the file: ", trim(file_name)
    call Cg_Error_Exit_F()
  endif

  end subroutine