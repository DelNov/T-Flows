!==============================================================================!
  subroutine File_Mod_Set_Name(name_out, time_step, processor, extension)
!------------------------------------------------------------------------------!
!   Creates the file name depending on time step, subdomain and file type.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: name_out
  integer, optional :: time_step
  integer, optional :: processor
  character(len=*)  :: extension
!-----------------------------------[Locals]-----------------------------------!
  integer          :: last_pos
  integer          :: lext
  character(len=5) :: num_proc  ! processor number as a string
!==============================================================================!

  name_out = problem_name
  last_pos = len_trim(name_out)

  !----------------------------------!
  !   Append time step to name_out   !
  !----------------------------------!
  if(present(time_step)) then
    write(name_out(last_pos+1:last_pos+3), '(a3)')   '-ts'
    write(name_out(last_pos+4:last_pos+9), '(i6.6)') time_step
    last_pos = last_pos + 9
  end if

  !-----------------------------!
  !   Append processor number   !
  !-----------------------------!
  if(present(processor)) then
    if(processor > 0) then
      write(name_out(last_pos+1:last_pos+3), '(a3)')   '-pu'
      write(name_out(last_pos+4:last_pos+8), '(i5.5)') processor
      last_pos = last_pos + 8
    end if
  end if

  !--------------------------!
  !   Append the extension   !
  !--------------------------!
  lext = len_trim(extension)
  name_out(last_pos+1:last_pos+lext) = extension(1:lext)

  end subroutine
