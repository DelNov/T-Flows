!==============================================================================!
  subroutine Grid_Mod_Sort_Cells_Smart(grid)
!------------------------------------------------------------------------------!
!   Sorts cells by their geometrical positions.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: s, c, c1, c2
  integer, allocatable :: old_nn   (:)    ! old number of nodes (per cell)
  integer, allocatable :: old_nodes(:,:)  ! old nodes list per cell
  integer, allocatable :: old_nc   (:)    ! old number of cells (per cell)
  integer, allocatable :: old_nf   (:)    ! old number of faces (per cell)
  integer, allocatable :: old_faces(:,:)  ! old faces list per cell
  integer, allocatable :: old_cells(:,:)  ! old cells list per cell
  real, allocatable    :: xc(:), yc(:), zc(:)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: L = MAX_CELLS_N_CELLS
  integer, parameter :: M = MAX_CELLS_N_NODES
  integer, parameter :: N = MAX_CELLS_N_FACES
!==============================================================================!

  allocate(old_nn      (grid % n_cells));  old_nn   (:)   = 0
  allocate(old_nodes(M, grid % n_cells));  old_nodes(:,:) = 0
  allocate(old_nc      (grid % n_cells));  old_nc   (:)   = 0
  allocate(old_cells(L, grid % n_cells));  old_cells(:,:) = 0
  allocate(old_nf      (grid % n_cells));  old_nf   (:)   = 0
  allocate(old_faces(N, grid % n_cells));  old_faces(:,:) = 0
  allocate(xc          (grid % n_cells));  xc       (:)   = 0.0
  allocate(yc          (grid % n_cells));  yc       (:)   = 0.0
  allocate(zc          (grid % n_cells));  zc       (:)   = 0.0

  !------------------------!
  !   Store cells' nodes   !
  !------------------------!
  do c = 1, grid % n_cells
    old_nn   (c)      = grid % cells_n_nodes(c)
    old_nodes(1:M, c) = grid % cells_n(1:M, c)
    old_nf   (c)      = grid % cells_n_faces(c)
    old_faces(1:N, c) = grid % cells_f(1:N, c)
    old_nc   (c)      = grid % cells_n_cells(c)
    old_cells(1:L, c) = grid % cells_c(1:L, c)
  end do

  !--------------------------!
  !   Set sorting criteria   !
  !--------------------------!
  do c = 1, grid % n_cells
    xc(c) = grid % xc(c)
    yc(c) = grid % yc(c)
    zc(c) = grid % zc(c)
    grid % old_c(c) = c
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort_Mod_3_Real_Carry_Int(xc(1:grid % n_cells),  &
                                 yc(1:grid % n_cells),  &
                                 zc(1:grid % n_cells),  &
                                 grid % old_c(1:grid % n_cells))
  ! This is a bit of a bluff
  do c = 1, grid % n_cells
    grid % new_c(grid % old_c(c)) = c
  end do

  !----------------------------------!
  !   Update cell numbers at faces   !
  !----------------------------------!
  do s = 1, grid % n_faces + grid % n_shadows
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)
    grid % faces_c(1, s) = grid % new_c(c1)
    if(c2 > 0) then
      grid % faces_c(2, s) = grid % new_c(c2)

      ! If the face changed its orientation during cell renumeration
      if(grid % faces_c(2, s) < grid % faces_c(1, s)) then
        call Swap_Int(grid % faces_c(1, s), grid % faces_c(2, s))
        grid % sx(s) = -grid % sx(s)
        grid % sy(s) = -grid % sy(s)
        grid % sz(s) = -grid % sz(s)
        grid % dx(s) = -grid % dx(s)
        grid % dy(s) = -grid % dy(s)
        grid % dz(s) = -grid % dz(s)
        grid % f (s) = 1.0 - grid % f(s)
      end if
    end if
  end do

  !-----------------------------------------------!
  !   Do the sorting of data pertinent to cells   !
  !-----------------------------------------------!
  do c = 1, grid % n_cells
    grid % cells_n_nodes(grid % new_c(c)) = old_nn   (c)
    grid % cells_n(1:M,  grid % new_c(c)) = old_nodes(1:M, c)
    grid % cells_n_faces(grid % new_c(c)) = old_nf   (c)
    grid % cells_f(1:N,  grid % new_c(c)) = old_faces(1:N, c)
    grid % cells_n_cells(grid % new_c(c)) = old_nc   (c)
    grid % cells_c(1:N,  grid % new_c(c)) = old_cells(1:N, c)
  end do
  call Sort_Mod_Real_By_Index(grid % n_cells, grid % xc (1), grid % new_c(1))
  call Sort_Mod_Real_By_Index(grid % n_cells, grid % yc (1), grid % new_c(1))
  call Sort_Mod_Real_By_Index(grid % n_cells, grid % zc (1), grid % new_c(1))
  call Sort_Mod_Real_By_Index(grid % n_cells, grid % vol(1), grid % new_c(1))

  end subroutine
