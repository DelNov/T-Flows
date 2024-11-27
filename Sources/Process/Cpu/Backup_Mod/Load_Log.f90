!==============================================================================!
  subroutine Load_Log(Backup, Comm, disp, vc, var_name, var_value)
!------------------------------------------------------------------------------!
!   Reads a single named logical variable from backup file.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type) :: Backup
  type(Comm_Type)    :: Comm
  integer(DP)        :: disp
  integer            :: vc
  character(len=*)   :: var_name
  logical            :: var_value
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: vn
  integer       :: vs, cnt_loop
  integer(DP)   :: disp_loop
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Backup)
!==============================================================================!

  cnt_loop  = 0
  disp_loop = 0

  !--------------------------------------------------------!
  !   Browse the entire file until you find the variable   !
  !--------------------------------------------------------!
  do

    ! Increase counter
    cnt_loop = cnt_loop + 1

    call Comm % Read_Text(fh, vn, disp_loop)  ! variable name
    call Comm % Read_Int (fh, vs, disp_loop)  ! variable offset

    ! If variable is found, read it and retrun
    if(vn .eq. var_name) then
      if(First_Proc()) print *, '# Reading variable: ', trim(vn)
      call Comm % Read_Log(fh, var_value, disp_loop)
      disp = disp_loop
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

    ! Check if variable is in the file
    if(cnt_loop >= vc) goto 1

  end do

1 if(First_Proc()) print *, '# Variable: ', trim(var_name), ' not found! ',  &
                             'Setting the value to false!'
  var_value = .false.

  end subroutine
