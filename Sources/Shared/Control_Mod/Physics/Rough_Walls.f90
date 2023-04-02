!==============================================================================!
  subroutine Rough_Walls(Control, rough, verbose)
!------------------------------------------------------------------------------!
!   Reading wall roughness from the control file.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  logical, intent(out) :: rough
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('ROUGH_WALLS', 'no', val, verbose)
  call String % To_Upper_Case(val)

  select case(val)

    case('YES')
      rough = .true.

    case('NO')
      rough = .false.

    case default
      call Message % Error(60,                                            &
                           'Unknown wall roughness state: '//trim(val)//  &
                           '. \n This error is critical. Exiting.',       &
                           file=__FILE__, line=__LINE__, one_proc=.true.)
  end select

  end subroutine
