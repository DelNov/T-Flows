!==============================================================================!
  subroutine Balanced_Pressure_Gradient(Control, bal, def, verbose)
!------------------------------------------------------------------------------!
!>  Reads from the control file if cell-centered pressure gradients are to be
!>  computed in balance with body forces (important for stratified flows).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  logical, intent(out) :: bal      !! output value, true or false
  logical, intent(in)  :: def      !! default value
  logical, optional    :: verbose  !! verbosity of the output
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val, def_val
!==============================================================================!

  def_val = 'no'
  if(def) def_val = 'yes'

  call Control % Read_Char_Item('BALANCED_PRESSURE_GRADIENT', def_val,  &
                                val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    bal = .true.

  else if( val .eq. 'NO' ) then
    bal = .false.

  else
    call Message % Error(60,                                              &
             'Unknown state for balanced pressure gradient: '//trim(val)//&
             '. \n This error is critical.  Exiting.',                    &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
