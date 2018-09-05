!==============================================================================!
  subroutine Control_Mod_Turbulence_Statistics(verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence statistics from the control file                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turbulence_Mod, only: turbulence_statistics
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENCE_STATISTICS', 'no',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if     (val .eq. 'YES') then
    turbulence_statistics = .true.

  else if(val .eq. 'NO') then
    turbulence_statistics = .false.

  else
    print *, '# Unknown optoin for turbuelnce statistics: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine
