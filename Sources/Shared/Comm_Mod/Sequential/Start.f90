!==============================================================================!
  subroutine Comm_Mod_Start
!------------------------------------------------------------------------------!
!   Initializes sequential execution.                                          !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer              :: n
  integer, allocatable :: seeds(:)
!==============================================================================!

  Communicator % this_processor = 0
  Communicator % n_processors   = 0

  call random_seed(size = n)
  allocate(seeds(n)); seeds(:) = 12321
  call random_seed(put = seeds)

  end subroutine
