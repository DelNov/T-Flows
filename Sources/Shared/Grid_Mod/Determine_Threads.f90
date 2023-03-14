!==============================================================================!
  subroutine Determine_Threads(Grid)
!------------------------------------------------------------------------------!
!   Take care of threads for OpenMP                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer :: c_f, c_l, s_f, s_l      ! first/last cell, first/last face
!==============================================================================!

  !----------------------------!
  !   Find number of threads   !
  !----------------------------!
  Grid % n_threads = 0
# ifdef _OPENMP
    Grid % n_threads = omp_get_max_threads()
# endif

  !----------------------------------------------------------!
  !   If there is more than one thread, decompose the mesh   !
  !----------------------------------------------------------!
  if(PROGRAM_NAME == "Process" .and.  &
     Grid % n_threads > 1) then

    ! Set the face range to the inside faces
    s_f = Grid % region % f_face(Grid % n_regions)
    s_l = Grid % region % l_face(Grid % n_regions)

    ! Create METIS
    call Metis % Create_Metis(s_f, s_l, Grid % faces_c, Grid % n_threads)

    ! Set the cell range to the inside faces
    c_f = Grid % region % f_cell(Grid % n_regions)
    c_l = Grid % region % l_cell(Grid % n_regions)

    ! Call METIS
    call Metis % Call_Metis(Grid % n_threads,  &
                            Grid % cell_thread(c_f:c_l))
  end if

  end subroutine
