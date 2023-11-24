!==============================================================================!
  subroutine Set_Num_Threads(Omp, n)
!------------------------------------------------------------------------------!
!>  The Set_Num_Threads subroutine in the Omp_Mod module serves to specify the
!> number of threads to be used for OpenMP parallel processing.
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Omp_Type)     :: Omp  !! parent class
  integer, intent(in) :: n    !! number of threads
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Omp)
!==============================================================================!

# ifdef _OPENMP
    call omp_set_num_threads(n)
# else
    Unused(n)
# endif

  end subroutine

