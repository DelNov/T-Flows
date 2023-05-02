!==============================================================================!
  logical function Time_To_Save_Backup(Backup)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type) :: Backup
!==============================================================================!

  Time_To_Save_Backup = mod(Time % Curr_Dt(), backup % interval) .eq. 0  &
                            .and. .not. Time % Curr_Dt() .eq. 0

  end function
