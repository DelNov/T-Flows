!==============================================================================!
  subroutine Name_File(sub, name_out, ext)
!------------------------------------------------------------------------------!
!   Creates the file name depending on the subdomain and file type.            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: sub
  character(len=*)  :: name_out
  character(len=*)  :: ext
!-----------------------------------[Locals]-----------------------------------!
  integer          :: lext
  character(len=5) :: num_proc  ! processor number as a string
!==============================================================================!

  ! Take the extension length
  lext = len_trim(ext)

  ! Set the stem to problem name
  name_out = problem_name

  ! For parallel runs, add processor number
  if(sub > 0) then
    write(num_proc,'(i5.5)') sub
    write(name_out(len_trim(name_out)+1:len_trim(name_out)+3),'(a3)') '-pu' 
    write(name_out(len_trim(name_out)+1:len_trim(name_out)+5),'(a5)') num_proc
  end if 

  ! Add file extension
  name_out(len_trim(name_out)+1:len_trim(name_out)+lext) = ext(1:lext) 

  end subroutine
