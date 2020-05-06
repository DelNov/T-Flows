!==============================================================================!
  subroutine File_Mod_Set_Name(name_out,   &
                               time_step,  &
                               processor,  &
                               appendix,   &
                               extension,  &
                               domain)
!------------------------------------------------------------------------------!
!   Creates the file name depending on time step, subdomain, type and domain   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)           :: name_out
  integer,          optional :: time_step
  integer,          optional :: processor
  character(len=*), optional :: appendix   ! used to add '-bnd' to name
  character(len=*)           :: extension
  integer,          optional :: domain
!-----------------------------------[Locals]-----------------------------------!
  integer          :: last_pos
  integer          :: lext, lapp
  character(len=5) :: num_proc  ! processor number as a string
!==============================================================================!

  !-------------------------------------------!
  !   Handle problems with multiple domains   !
  !-------------------------------------------!
  if(present(domain)) then
    name_out = problem_name(domain)
  else
    name_out = problem_name(1)
  end if

  last_pos = len_trim(name_out)

  !----------------------------------!
  !   Add appendix to problem name   !
  !----------------------------------!
  if(present(appendix)) then
    lapp = len_trim(appendix)
    name_out(last_pos+1:last_pos+lapp) = appendix(1:lapp)
    last_pos = last_pos + lapp
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

  !----------------------------------!
  !   Append time step to name_out   !
  !----------------------------------!
  if(present(time_step)) then
    write(name_out(last_pos+1:last_pos+3), '(a3)')   '-ts'
    write(name_out(last_pos+4:last_pos+9), '(i6.6)') time_step
    last_pos = last_pos + 9
  end if

  !--------------------------!
  !   Append the extension   !
  !--------------------------!
  lext = len_trim(extension)
  name_out(last_pos+1:last_pos+lext) = extension(1:lext)

  end subroutine
