!==============================================================================!
  subroutine Start_Parallel(Global)
!------------------------------------------------------------------------------!
!   Initializes sequential execution.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Global
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
