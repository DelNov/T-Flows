!==============================================================================!
  subroutine Sgs_Buoyancy_Correction(Control, corr, verbose)
!------------------------------------------------------------------------------!
!>  Reads if the turbulent viscosity of LES models is reduced in stable
!>  stratification, with sqrt(1 - Ri / Pr_t) after Eidson (JFM, 1985).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  logical, intent(out) :: corr     !! true if correction is engaged
  logical, optional    :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('SGS_BUOYANCY_CORRECTION', 'yes', val, verbose)
  call String % To_Upper_Case(val)

  select case(val)

    case('YES')
      corr = .true.

    case('NO')
      corr = .false.

    case default
      call Message % Error(60,                                              &
                           'Unknown state of SGS buoyancy correction: '//   &
                           trim(val)//'. \n This error is critical. Exiting.',&
                           file=__FILE__, line=__LINE__, one_proc=.true.)
  end select

  end subroutine
