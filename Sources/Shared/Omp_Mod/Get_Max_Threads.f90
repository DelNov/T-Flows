!==============================================================================!
  integer function Get_Max_Threads(Omp)
!------------------------------------------------------------------------------!
!>  The Get_Max_Threads function in the Omp_Mod module is designed to obtain
!>  the maximum number of threads available for OpenMP parallel processing.
!>  If the code is compiled without OpenMP, it will select only one thread.
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Omp_Type) :: Omp  !! parent class
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Omp)
!==============================================================================!

  Get_Max_Threads = 1
# ifdef _OPENMP
    Get_Max_Threads = omp_get_max_threads()
# endif

  end function

