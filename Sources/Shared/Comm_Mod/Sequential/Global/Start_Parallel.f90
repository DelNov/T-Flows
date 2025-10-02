!==============================================================================!
  subroutine Start_Parallel(Global)
!------------------------------------------------------------------------------!
!>  Initializes a sequential execution context. This subroutine sets up the
!>  necessary environment for running a program in a non-parallel (single
!>  processor) mode. It configures the global communication object for a
!>  sequential run and initializes the random number generator with a fixed
!>  seed to ensure reproducibility.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(out) :: Global  !! global communication object
!-----------------------------------[Locals]-----------------------------------!
  integer              :: n
  integer, allocatable :: seeds(:)
!==============================================================================!

  Global % this_processor = 0
  Global % n_processors   = 0

  call random_seed(size = n)
  allocate(seeds(n)); seeds(:) = 12321
  call random_seed(put = seeds)

  end subroutine
