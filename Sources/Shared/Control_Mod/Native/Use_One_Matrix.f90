!==============================================================================!
  subroutine Use_One_Matrix(Control, use_one, verbose)
!------------------------------------------------------------------------------!
!>  Reads skewness correction for VOF.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control    !! parent class
  logical, intent(out) :: use_one
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('USE_ONE_MATRIX', 'yes', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    use_one = .true.

  else if( val .eq. 'NO' ) then
    use_one = .false.

  else
    call Message % Error(72,                                    &
             'Unknown state for use one matrix: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',          &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
