!==============================================================================!
  subroutine End_Parallel(Global)
!------------------------------------------------------------------------------!
!>  Terminates parallel execution in a distributed computing environment.
!>  This subroutine is crucial for closing parallel communication and ensuring
!>  a synchronized exit of all processes in a multi-processor setting.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Global  !! global communication object
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  call Mpi_Barrier(MPI_COMM_WORLD, error)
  call Mpi_Finalize(error)

  end subroutine
