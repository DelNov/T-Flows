!==============================================================================!
  subroutine Start_Parallel(Global)
!------------------------------------------------------------------------------!
!>  Initializes the parallel execution environment in a distributed computing
!>  system using MPI (Message Passing Interface). This subroutine sets up the
!>  essential components for parallel processing, including the total number
!>  of processors, the rank of the current processor, and initializes random
!>  seeds for consistent random number generation across different processors.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(out) :: Global  !! global communication object
!-----------------------------------[Locals]-----------------------------------!
  integer              :: error
  integer              :: n
  integer, allocatable :: seeds(:)
!==============================================================================!

  call Mpi_Init(error)

  ! Get number of processors
  call Mpi_Comm_Size(MPI_COMM_WORLD, Global % n_processors, error)

  ! Get current processor number
  call Mpi_Comm_Rank(MPI_COMM_WORLD, Global % this_processor, error)

  ! Use Fortran counting - from 1
  Global % this_processor = Global % this_processor + 1

  ! But if run is sequential, set the only processor to zero
  if(Global % n_processors .eq. 1) then
    Global % n_processors   = 0
    Global % this_processor = 0
  endif

  ! Set proper types for communication
  comm_type_int  = MPI_INTEGER
  comm_type_log  = MPI_LOGICAL
  comm_type_real = MPI_DOUBLE_PRECISION
  if(RP .eq. SP) then
    comm_type_real = MPI_REAL
  end if

  ! Make sure that integers are 32-bit
  if(IP .eq. DP) then
    if(Global % this_processor < 2) then
      print *, '# Error - 64 bit integers are not supported!'
      print *, '# This error is critical, exiting!'
      call Global % End_Parallel
      stop
    end if
  end if

  ! Initialize random seeds (to make all processors
  ! generate the same sequence of random numbers)
  call random_seed(size = n)
  allocate(seeds(n)); seeds(:) = 12321
  call random_seed(put = seeds)

  end subroutine
