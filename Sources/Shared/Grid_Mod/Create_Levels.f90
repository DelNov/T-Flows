!==============================================================================!
  subroutine Grid_Mod_Create_Levels(grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: lev, n_cells, n_faces, c, s
!==============================================================================!

  !--------------------------------!
  !                                !
  !   Allocate memory for levels   !
  !                                !
  !--------------------------------!

  ! At each level, allocate memory for individual data
  do lev = grid % n_levels + 1, 1, -1

    ! Tentative number of cells and faces
    n_cells = N_MG_PARTS ** (grid % n_levels + 1 - lev)
    n_faces = grid % n_faces

    if(lev .eq. 1) then
      n_cells = grid % n_cells
      n_faces = grid % n_faces
    end if

    ! Connectivity with finest grid
    allocate(grid % level(lev) % cell(grid % n_cells))
    allocate(grid % level(lev) % face(grid % n_faces))

    ! Intra-level connectivity and geometry
    allocate(grid % level(lev) % faces_c(2, n_faces))
    allocate(grid % level(lev) % n_finest_cells(n_cells))
    allocate(grid % level(lev) % xc     (   n_cells))
    allocate(grid % level(lev) % yc     (   n_cells))
    allocate(grid % level(lev) % zc     (   n_cells))
  end do

  !---------------------------------------------------------!
  !                                                         !
  !   Initialize arrays for level zero and coarser levels   !
  !                                                         !
  !---------------------------------------------------------!

  !------------------------------------!
  !     Cells and faces on level 1     !
  !   (Level 1 is the original grid)   !
  !------------------------------------!

  ! Cells
  grid % level(1) % n_cells = grid % n_cells
  do c = 1, grid % n_cells
    grid % level(1) % cell(c) = c
  end do

  ! Faces
  grid % level(1) % n_faces = grid % n_faces
  do s = 1, grid % n_faces
    grid % level(1) % face(s) = s
    grid % level(1) % faces_c(1, s) = grid % faces_c(1, s)
    grid % level(1) % faces_c(2, s) = grid % faces_c(2, s)
  end do

  !-------------------------------------!
  !   Cells and faces on other levels   !
  !-------------------------------------!
  do lev = 2, grid % n_levels + 1

    ! Cells
    do c = 1, grid % n_cells
      grid % level(lev) % cell(c) = 1
    end do

    ! Faces
    do s = 1, grid % n_faces
      grid % level(lev) % face(s) = 0
    end do
  end do

  end subroutine
