!==============================================================================!
  subroutine Comm_Mod_Start
!------------------------------------------------------------------------------!
!   Initializes parallel execution.                                            !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  call Mpi_Init(error)

  ! Get number of processors
  call Mpi_Comm_Size(MPI_COMM_WORLD, n_proc, error) 

  ! Get current processor number
  call Mpi_Comm_Rank(MPI_COMM_WORLD, this_proc, error)

  ! Use Fortran counting - from 1
  this_proc = this_proc + 1

  ! But if run is sequential, set the only processor to zero
  if(n_proc .eq. 1) then
    n_proc = 0
    this_proc = 0
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
    if(this_proc < 2) then
      print *, '# Error - 64 bit integers are not supported!'
      print *, '# This error is critical, exiting!'
      call Comm_Mod_End
      stop
    end if
  end if

  end subroutine
