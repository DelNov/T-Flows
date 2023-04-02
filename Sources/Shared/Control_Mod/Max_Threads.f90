!==============================================================================!
  subroutine Max_Threads(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

# ifdef _OPENMP
  call Control % Read_Int_Item('MAX_THREADS', 2, val, verbose)
# else
  val = 1
  Unused(verbose)
# endif

  end subroutine
