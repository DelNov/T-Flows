!==============================================================================!
  subroutine Monin_Obukov(Control, most, verbose)
!------------------------------------------------------------------------------!
!>  Reads the Monin-Obukov similarity theory (MOST) model from control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  logical, intent(out) :: most     !! true if MOST model is to be used
  logical, optional    :: verbose  !! controls output verbosity
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
      call Message % Error(60,                                       &
                           'Unknown ABL model: '//trim(val)//        &
                           '. \n This error is critical. Exiting.',  &
                           file=__FILE__, line=__LINE__, one_proc=.true.)
  end select

  end subroutine
