!==============================================================================!
  subroutine Set_Name(File,       &
                      name_out,   &
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
  class(File_Type)           :: File
  character(len=*)           :: name_out
  integer,          optional :: time_step
  integer,          optional :: processor
  character(len=*), optional :: appendix   ! used to add '-bnd' to name
  character(len=*)           :: extension
  integer,          optional :: domain
!-----------------------------------[Locals]-----------------------------------!
  integer       :: last_pos
  integer       :: ldir, lnam, lext, lapp
  character(SL) :: rel_path, sys_comm
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  !-------------------------------------------!
  !   Create directories for each processor   !
  !-------------------------------------------!
  last_pos = 0
  if(present(processor)) then
    if(processor > 0) then
      rel_path = 'Sub/00000/'
      ldir = len_trim(rel_path)
      write(rel_path(ldir-5:ldir-1), '(i5.5)') processor
      sys_comm = 'mkdir -p ' // trim(rel_path)
      call system(trim(sys_comm))
      name_out = rel_path
      last_pos = len_trim(name_out)
    end if
  end if

  !-------------------------------------------!
  !   Handle problems with multiple domains   !
  !-------------------------------------------!
  if(present(domain)) then
    Assert(domain >=  1)
    Assert(domain <= MD)
    lnam = len_trim(problem_name(domain))
    if(last_pos > 0) then
      name_out(last_pos+1:last_pos+lnam) = problem_name(domain)(1:lnam)
    else
      name_out = problem_name(domain)
    end if
  else
    lnam = len_trim(problem_name(1))
    if(last_pos > 0) then
      name_out(last_pos+1:last_pos+lnam) = problem_name(1)(1:lnam)
    else
      name_out = problem_name(1)
    end if
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
