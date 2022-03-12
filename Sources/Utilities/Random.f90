!==============================================================================!
  program Random
!------------------------------------------------------------------------------!
!   Checks how random_number and random_seed work in parallel                  !
!                                                                              !
!   Compile with:  mpif90 -o Random Random.f90                                 !
!                                                                              !
!   Exectute as:   mpirun -np 4 ./Random                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Mpi_f08
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: error, n_proc, this_proc, n, i
  integer, allocatable :: seeds(:)
  real    :: numbah
!==============================================================================!

  ! Initialize MPI
  call Mpi_Init(error)

  ! Get number of processors
  call Mpi_Comm_Size(MPI_COMM_WORLD, n_proc, error) 

  ! Get current processor number
  call Mpi_Comm_Rank(MPI_COMM_WORLD, this_proc, error)

  ! Use Fortran counting - from 1
  this_proc = this_proc + 1

  call random_seed(size = n)
  allocate(seeds(n)); seeds(:) = 0
  call random_seed(put = seeds)

  do i = 1, 20
    call random_number(numbah)
    print *, '# Random number from', this_proc, 'is:', numbah
  end do

  ! End MPI
  call Mpi_Finalize(error)

  end program
