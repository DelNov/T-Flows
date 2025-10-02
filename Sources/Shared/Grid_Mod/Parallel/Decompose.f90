!==============================================================================!
  subroutine Decompose(Grid, n_parts)
!------------------------------------------------------------------------------!
!>  Decomposes the grid with METIS library.  Cleraly, it is called from Divide.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid     !! grid to be decomposed
  integer,          intent(in)    :: n_parts  !! number of sub-domains
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, s, s_f, s_l
  integer, allocatable :: part(:)          ! result of partitioning
!==============================================================================!

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  allocate(part(Grid % n_cells))

  !----------------------------------------------------------------------------!
  !   Prepare METIS and Exectute a call to METIS graph partitioning function   !
  !----------------------------------------------------------------------------!

  ! Set the bounding faces
  s_l = Grid % n_faces + Grid % n_shadows  ! last face at this point
  do s_f = 1, s_l                          ! search for first inside face
    if(Grid % faces_c(2, s_f) > 0) exit    ! when you find first inside face
  end do

  ! Call METIS with estimated face range
  call Metis % Create_Metis(s_f, s_l, Grid % faces_c, n_parts)

  !-----------------------------!
  !   Execute a call to METIS   !
  !-----------------------------!
  call Metis % Call_Metis(n_parts, part)

  !-----------------------------------------------------!
  !   Save the result from the call to METIS function   !
  !-----------------------------------------------------!
  do c = 1, Grid % n_cells
    Grid % Comm % cell_proc(c) = part(c)
  end do

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 < 0) then
      Grid % Comm % cell_proc(c2) = Grid % Comm % cell_proc(c1)
    end if
  end do

  end subroutine
