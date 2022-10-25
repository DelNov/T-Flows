!==============================================================================!
  subroutine Sort_Cells_Smart(Grid)
!------------------------------------------------------------------------------!
!   Sorts cells by their geometrical positions.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
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
  integer :: mcc, mcn, mcf
!==============================================================================!

  mcc = size(Grid % cells_c, 1)
  mcn = size(Grid % cells_n, 1)
  mcf = size(Grid % cells_f, 1)

  allocate(old_nn        (Grid % n_cells));  old_nn   (:)   = 0
  allocate(old_nodes(mcn, Grid % n_cells));  old_nodes(:,:) = 0
  allocate(old_nc        (Grid % n_cells));  old_nc   (:)   = 0
  allocate(old_cells(mcc ,Grid % n_cells));  old_cells(:,:) = 0
  allocate(old_nf        (Grid % n_cells));  old_nf   (:)   = 0
  allocate(old_faces(mcf, Grid % n_cells));  old_faces(:,:) = 0
  allocate(xc            (Grid % n_cells));  xc       (:)   = 0.0
  allocate(yc            (Grid % n_cells));  yc       (:)   = 0.0
  allocate(zc            (Grid % n_cells));  zc       (:)   = 0.0

  !------------------------!
  !   Store cells' nodes   !
  !------------------------!
  do c = 1, Grid % n_cells
    old_nn   (c)        = Grid % cells_n_nodes(c)
    old_nodes(1:mcn, c) = Grid % cells_n(1:mcn, c)
    old_nf   (c)        = Grid % cells_n_faces(c)
    old_faces(1:mcf, c) = Grid % cells_f(1:mcf, c)
    old_nc   (c)        = Grid % cells_n_cells(c)
    old_cells(1:mcc, c) = Grid % cells_c(1:mcc, c)
  end do

  !--------------------------!
  !   Set sorting criteria   !
  !--------------------------!
  do c = 1, Grid % n_cells
    xc(c) = Grid % xc(c)
    yc(c) = Grid % yc(c)
    zc(c) = Grid % zc(c)
    Grid % old_c(c) = c
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort % Three_Real_Carry_Int(xc(1:Grid % n_cells),  &
                                   yc(1:Grid % n_cells),  &
                                   zc(1:Grid % n_cells),  &
                                   Grid % old_c(1:Grid % n_cells))
  ! This is a bit of a bluff
  do c = 1, Grid % n_cells
    Grid % new_c(Grid % old_c(c)) = c
  end do

  !----------------------------------!
  !   Update cell numbers at faces   !
  !----------------------------------!
  do s = 1, Grid % n_faces + Grid % n_shadows
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    Grid % faces_c(1, s) = Grid % new_c(c1)
    if(c2 > 0) then
      Grid % faces_c(2, s) = Grid % new_c(c2)

      ! If the face changed its orientation during cell renumeration
      if(Grid % faces_c(2, s) < Grid % faces_c(1, s)) then
        call Swap_Int(Grid % faces_c(1, s), Grid % faces_c(2, s))
        Grid % sx(s) = -Grid % sx(s)
        Grid % sy(s) = -Grid % sy(s)
        Grid % sz(s) = -Grid % sz(s)
        Grid % dx(s) = -Grid % dx(s)
        Grid % dy(s) = -Grid % dy(s)
        Grid % dz(s) = -Grid % dz(s)
        Grid % f (s) = 1.0 - Grid % f(s)
      end if
    end if
  end do

  !-----------------------------------------------!
  !   Do the sorting of data pertinent to cells   !
  !-----------------------------------------------!
  do c = 1, Grid % n_cells
    Grid % cells_n_nodes(Grid % new_c(c)) = old_nn          (c)
    Grid % cells_n(1:mcn,Grid % new_c(c)) = old_nodes(1:mcn, c)
    Grid % cells_n_faces(Grid % new_c(c)) = old_nf          (c)
    Grid % cells_f(1:mcf,Grid % new_c(c)) = old_faces(1:mcf, c)
    Grid % cells_n_cells(Grid % new_c(c)) = old_nc          (c)
    Grid % cells_c(1:mcc,Grid % new_c(c)) = old_cells(1:mcc, c)
  end do
  call Sort % Real_By_Index(Grid % n_cells, Grid % xc (1), Grid % new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % yc (1), Grid % new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % zc (1), Grid % new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % vol(1), Grid % new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % ixx(1), Grid % new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % iyy(1), Grid % new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % izz(1), Grid % new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % ixy(1), Grid % new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % ixz(1), Grid % new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % iyz(1), Grid % new_c(1))

  end subroutine
