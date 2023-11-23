!==============================================================================!
  subroutine Max_Threads(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads number of threads to use (for OMP) from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! number of threads to use
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

# ifdef _OPENMP
  call Control % Read_Int_Item('MAX_THREADS', 2, val, verbose)
# else
  val = 1
  Unused(Control)
  Unused(verbose)
# endif

  end subroutine
