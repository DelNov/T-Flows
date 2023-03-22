!==============================================================================!
  logical function Time_To_Save_Backup(Backup, curr_dt)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type) :: Backup
  integer            :: curr_dt  ! current time step
!==============================================================================!

  Time_To_Save_Backup = mod(curr_dt, backup % interval) .eq. 0  &
                            .and. .not. curr_dt .eq. 0

  end function
