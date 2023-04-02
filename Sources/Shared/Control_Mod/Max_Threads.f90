!==============================================================================!
  subroutine Control_Mod_Max_Threads(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
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
