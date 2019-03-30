!==============================================================================!
  subroutine Backup_Mod_Read_Int_Array(fh, disp, vc, arr_name, arr_value)
!------------------------------------------------------------------------------!
!   Reads a named integer array from backup file.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer               :: fh, disp, vc
  character(len=*)      :: arr_name
  integer, dimension(:) :: arr_value
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: vn
  integer           :: vs, disp_loop, cnt_loop, length
!==============================================================================!

  cnt_loop  = 0
  disp_loop = 0
  length    = size(arr_value)

  !--------------------------------------------------------!
  !   Browse the entire file until you find the variable   !
  !--------------------------------------------------------!
  do

    ! Increase counter
    cnt_loop = cnt_loop + 1

    call Comm_Mod_Read_Text(fh, vn, disp_loop)  ! variable name
    call Comm_Mod_Read_Int (fh, vs, disp_loop)  ! variable offset

    ! If variable is found, read it and retrun
    if(vn .eq. arr_name) then
      if(this_proc < 2) print *, '# Reading array: ', trim(vn)
      call Comm_Mod_Read_Int_Array(fh, arr_value, disp_loop)
      disp = disp_loop
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

    ! Check if variable is in the file
    if(cnt_loop >= vc) goto 1

  end do

1 if(this_proc < 2) print *, '# Array: ', trim(arr_name), ' not found! ',  &
                             'Setting the values to 0!'
  arr_value(:) = 0

  end subroutine
