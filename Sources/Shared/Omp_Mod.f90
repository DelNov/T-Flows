#include "../Shared/Unused.h90"

!==============================================================================!
  module Omp_Mod
!------------------------------------------------------------------------------!
!>  This module, Omp_Mod, is dedicated to OpenMP functionalities within the
!>  T-Flows project. It encapsulates common OpenMP operations, providing a
!>  structured and efficient approach to parallel processing on manycores.
!>  The module interfaces with Omp_Lib to access OpenMP functions and manages
!>  parallel regions through the Region_Mod.
!------------------------------------------------------------------------------!
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

  !--------------!
  !   OMP type   !
  !--------------!
  !> The Omp_Type defined within the Omp_Mod module is created to manage
  !> and facilitate OpenMP functionalities within the T-Flows program.
  type Omp_Type    ! used inside the Grid_Type

    ! Number of threads
    integer :: n_threads  !! actual number of OMP threads

    ! Desired number of threads read from control file
    integer :: d_threads  !! desired number of OMP threads

    type(Region_Type) :: thread  !! stores parallel regions in the code

    ! Thread for each cell and face
    integer, allocatable :: cell_thread(:)  !! map cells to threads
    integer, allocatable :: face_thread(:)  !! map faces to threads

    contains
      procedure :: Get_Max_Threads
      procedure :: Set_Num_Threads

  end type

  contains

#   include "Omp_Mod/Get_Max_Threads.f90"
#   include "Omp_Mod/Set_Num_Threads.f90"

  end module
