!==============================================================================!
  subroutine Monin_Obukov(Control, most, verbose)
!------------------------------------------------------------------------------!
!   Reading Monin-Obukov Similarity Theory (MOST) model from the control file. !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  logical, intent(out) :: most
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('MONIN_OBUKOV', 'no', val, verbose)
  call String % To_Upper_Case(val)

  select case(val)

    case('YES')
      most = .true.

    case('NO')
      most = .false.

    case default
      call Message % Error(60,                                            &
                           'Unknown ABL model: '//trim(val)//  &
                           '. \n This error is critical. Exiting.',       &
                           file=__FILE__, line=__LINE__, one_proc=.true.)
  end select

  end subroutine
