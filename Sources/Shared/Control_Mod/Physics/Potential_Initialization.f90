!==============================================================================!
  subroutine Control_Mod_Potential_Initialization(pot_init, verbose)
!------------------------------------------------------------------------------!
!   Reading potential initialization from the control file.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: pot_init
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('POTENTIAL_INITIALIZATION', 'no', val, verbose)
  call String % To_Upper_Case(val)

  select case(val)

    case('YES')
      pot_init = .true.

    case('NO')
      pot_init = .false.

    case default
      call Message % Error(72,                                              &
               'Unknown state for potential initialization: '//trim(val)//  &
               '. \n This error is critical.  Exiting.',                    &
               file=__FILE__, line=__LINE__, one_proc=.true.)
  end select

  end subroutine
