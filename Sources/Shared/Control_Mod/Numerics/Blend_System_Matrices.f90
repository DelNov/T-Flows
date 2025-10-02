!==============================================================================!
  subroutine Blend_System_Matrices(Control, blend, verbose)
!------------------------------------------------------------------------------!
!>  Reads if system matrices will be blended with advection.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  logical, intent(out) :: blend    !! output value, true or false
  logical,    optional :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('BLEND_SYSTEM_MATRICES', 'yes', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    blend = .true.

  else if( val .eq. 'NO' ) then
    blend = .false.

  else
    call Message % Error(60,                                           &
             'Unknown state for blend system matrices: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                 &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
