!==============================================================================!
  subroutine Save_Initial_Condition(Control, save_init_cond, verbose)
!------------------------------------------------------------------------------!
!>  Reads, from the control file, if you should save initial condition.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control         !! parent class
  logical             :: save_init_cond  !! should you save initial condition
  logical,   optional :: verbose         !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('SAVE_INITIAL_CONDITION', 'yes', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    save_init_cond = .true.

  else if( val .eq. 'NO' ) then
    save_init_cond = .false.

  else
    call Message % Error(72,                                            &
             'Unknown state for save initial condition: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                  &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
