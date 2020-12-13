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
  integer, allocatable :: new_c(:)
  integer, allocatable :: old_c(:)
  integer, allocatable :: old_nn   (:)    ! old number of nodes (per cell)
  integer, allocatable :: old_nodes(:,:)  ! old nodes list per cell
  integer, allocatable :: old_nf   (:)    ! old number of faces (per cell)
  integer, allocatable :: old_faces(:,:)  ! old faces list per cell
  real,    allocatable :: criteria (:,:)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: M = MAX_CELLS_N_NODES
  integer, parameter :: N = MAX_CELLS_N_FACES
!==============================================================================!

  allocate(old_nn   (   grid % n_cells));     old_nn   (:)   = 0
  allocate(old_nodes(M, grid % n_cells));     old_nodes(:,:) = 0
  allocate(old_nf   (   grid % n_cells));     old_nf   (:)   = 0
  allocate(old_faces(N, grid % n_cells));     old_faces(:,:) = 0
  allocate(new_c    (   grid % n_cells));     new_c(:)       = 0
  allocate(old_c    (   grid % n_cells));     old_c(:)       = 0
  allocate(criteria (   grid % n_cells, 3));  criteria(:,:)  = 0.0

  !------------------------!
  !   Store cells' nodes   !
  !------------------------!
  do c = 1, grid % n_cells
    old_nn   (c)      = grid % cells_n_nodes(c)
    old_nodes(1:M, c) = grid % cells_n(1:M, c)
    old_nf   (c)      = grid % cells_n_faces(c)
    old_faces(1:N, c) = grid % cells_f(1:N, c)
  end do

  !--------------------------!
  !   Set sorting criteria   !
  !--------------------------!
  do c = 1, grid % n_cells
    criteria(c, 1) = grid % xc(c)
    criteria(c, 2) = grid % yc(c)
    criteria(c, 3) = grid % zc(c)
    old_c(c)       = c
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort_Mod_3_Real_Carry_Int(criteria(1:grid % n_cells, 1),  &
                                 criteria(1:grid % n_cells, 2),  &
                                 criteria(1:grid % n_cells, 3),  &
                                 old_c   (1:grid % n_cells))
  ! This is a bit of a bluff
  do c = 1, grid % n_cells
    new_c(old_c(c)) = c
  end do

  !----------------------------------!
  !   Update cell numbers at faces   !
  !----------------------------------!
  do s = 1, grid % n_faces + grid % n_shadows
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)
    grid % faces_c(1, s) = new_c(c1)
    if(c2 > 0) then
      grid % faces_c(2, s) = new_c(c2)

      ! If the face changed its orientation during cell renumeration
      if(grid % faces_c(2, s) < grid % faces_c(1, s)) then
        call Swap_Int(grid % faces_c(1, s), grid % faces_c(2, s))
        grid % sx(s) = -grid % sx(s)
        grid % sy(s) = -grid % sy(s)
        grid % sz(s) = -grid % sz(s)
        grid % dx(s) = -grid % dx(s)
        grid % dy(s) = -grid % dy(s)
        grid % dz(s) = -grid % dz(s)
        grid % f (s) = 1.0 - grid % f (s)
      end if
    end if
  end do

  !-----------------------------------------------!
  !   Do the sorting of data pertinent to cells   !
  !-----------------------------------------------!
  do c = 1, grid % n_cells
    grid % cells_n_nodes(new_c(c)) = old_nn   (c)
    grid % cells_n(1:M, new_c(c))  = old_nodes(1:M, c)
    grid % cells_n_faces(new_c(c)) = old_nf   (c)
    grid % cells_f(1:N, new_c(c))  = old_faces(1:N, c)
  end do
  call Sort_Mod_Real_By_Index(grid % n_cells, grid % xc (1), new_c(1))
  call Sort_Mod_Real_By_Index(grid % n_cells, grid % yc (1), new_c(1))
  call Sort_Mod_Real_By_Index(grid % n_cells, grid % zc (1), new_c(1))
  call Sort_Mod_Real_By_Index(grid % n_cells, grid % vol(1), new_c(1))

  end subroutine
