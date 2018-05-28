!==============================================================================!
  subroutine Control_Mod_Pressure_Momentum_Coupling(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('PRESSURE_MOMENTUM_COUPLING', 'simple',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if(val .ne. 'SIMPLE'.and.  &
     val .ne. 'PROJECTION') then

    print *, '# Unknown pressure-momentum coupling: ', trim(val)
    print *, '# Exiting!'
    stop 

  end if

  end subroutine
