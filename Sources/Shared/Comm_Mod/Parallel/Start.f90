!==============================================================================!
  subroutine Comm_Mod_Start
!------------------------------------------------------------------------------!
!   Initializes parallel execution.                                            !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer              :: error
  integer              :: n
  integer, allocatable :: seeds(:)
!==============================================================================!

  call Mpi_Init(error)

  ! Get number of processors
  call Mpi_Comm_Size(MPI_COMM_WORLD, Communicator % n_processors, error)

  ! Get current processor number
  call Mpi_Comm_Rank(MPI_COMM_WORLD, Communicator % this_processor, error)

  ! Use Fortran counting - from 1
  Communicator % this_processor = Communicator % this_processor + 1

  ! But if run is sequential, set the only processor to zero
  if(Communicator % n_processors .eq. 1) then
    Communicator % n_processors   = 0
    Communicator % this_processor = 0
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
    if(Communicator % this_processor < 2) then
      print *, '# Error - 64 bit integers are not supported!'
      print *, '# This error is critical, exiting!'
      call Comm_Mod_End
      stop
    end if
  end if

  ! Initialize random seeds (to make all processors
  ! generate the same sequence of random numbers)
  call random_seed(size = n)
  allocate(seeds(n)); seeds(:) = 12321
  call random_seed(put = seeds)

  end subroutine
