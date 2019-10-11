!==============================================================================!
  subroutine Cgns_Mod_Open_File(file, mode)
!------------------------------------------------------------------------------!
!   Opens name_in file and set its file_id [parallel vesion]                   !
!------------------------------------------------------------------------------!
!   CGP_INDEPENDENT = any processor  CAN  access data                          !
!   CGP_COLLECTIVE  = all processors MUST access data                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*) :: file
  integer :: mode ! CG_MODE_READ, CG_MODE_WRITE or CG_MODE_MODIFY
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  file_name = file

  if (this_proc .lt. 2) print *, '# Opening file:', trim(file_name)

  ! Set the parallel IO mode for CGNS
  call Cgp_Pio_Mode_F(CGP_INDEPENDENT,  &
                      error)

  ! Open a CGNS file
  call Cgp_Open_F(file_name,  & !(in )
                  mode,       & !(in )
                  file_id,    & !(out)
                  error)        !(out)

  if (error .ne. 0) then
    print *, '# Failed to open the file: ', trim(file_name)
    call Cgp_Error_Exit_F()
  endif

  end subroutine
