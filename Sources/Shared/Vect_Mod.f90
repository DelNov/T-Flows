!==============================================================================!
  module Vect_Mod
!------------------------------------------------------------------------------!
!   Module for OpenMP functionality.                                           !
!                                                                              !
!   Some OpenMP functions you get with Omp_Lib:                                !
!                                                                              !
!   - int omp_get_max_threads()                                                !
!     Returns max possible (generally set by OMP_NUM_THREADS).                 !
!                                                                              !
!   - int omp_get_num_threads()                                                !
!     Returns number of threads in current team.                               !
!                                                                              !
!   - int omp_get_thread_num()                                                 !
!     Gets current thread number. It is between 0 and omp_get_num_threads()-1  !
!                                                                              !
!   - int omp_set_num_threads(num_threads)                                     !
!     Sets the number of threads overriding the OMP_NUM_THREADS.               !
!                                                                              !
!   I will try to encapsulate those as it is done in Comm_Mod/Type             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Omp_Lib       ! for OpenMP functionality
  use Region_Mod
  use Assert_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Vect type   !
  !---------------!
  type Vect_Type    ! used inside the Grid_Type

    ! Number of threads
    integer :: n_threads

    ! Desired number of threads read from control file
    integer :: d_threads

    type(Region_Type) :: region

    ! Thread i.d. defined for each cell
    integer, allocatable :: cell_thread(:)

    contains
      procedure :: Get_Max_Threads
      procedure :: Set_Num_Threads

  end type

  contains
  !==============================================================================!
    integer function Get_Max_Threads(Vect)
  !------------------------------------------------------------------------------!
    implicit none
  !------------------------------------------------------------------------------!
    class(Vect_Type) :: Vect
  !==============================================================================!

    Get_Max_Threads = 1
#   ifdef _OPENMP
      Get_Max_Threads = omp_get_max_threads()
#   endif

    end function

  !==============================================================================!
    subroutine Set_Num_Threads(Vect, n)
  !------------------------------------------------------------------------------!
    implicit none
  !------------------------------------------------------------------------------!
    class(Vect_Type) :: Vect
    integer, intent(in) :: n
  !==============================================================================!

#   ifdef _OPENMP
      call omp_set_num_threads(n)
#   endif

    end subroutine


  end module
