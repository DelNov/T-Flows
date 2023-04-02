!==============================================================================!
  subroutine Control_Mod_Rough_Walls(rough_walls, verbose)
!------------------------------------------------------------------------------!
!   Reading wall roughness from the control file.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: rough_walls
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('ROUGH_WALLS', 'no', val, verbose)
  call String % To_Upper_Case(val)

  select case(val)

    case('YES')
      rough_walls = .true.

    case('NO')
      rough_walls = .false.

    case default
      call Message % Error(60,                                            &
                           'Unknown wall roughness state: '//trim(val)//  &
                           '. \n This error is critical. Exiting.',       &
                           file=__FILE__, line=__LINE__, one_proc=.true.)
  end select

  end subroutine
