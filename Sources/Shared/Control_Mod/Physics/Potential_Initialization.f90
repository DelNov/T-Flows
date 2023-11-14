!==============================================================================!
  subroutine Potential_Initialization(Control, pot_init, verbose)
!------------------------------------------------------------------------------!
!>  Reads, from the control file, if potential initialization will be used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control   !! parent class
  logical, intent(out) :: pot_init  !! true if potential initialization is used
  logical, optional    :: verbose   !! controls output verbosity
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
