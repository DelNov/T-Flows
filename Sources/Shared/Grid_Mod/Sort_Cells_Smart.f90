!==============================================================================!
  Subroutine Grid_Mod_Sort_Cells_Smart(grid)
!------------------------------------------------------------------------------!
!   sorts array of faces in a smart way.  that would mean boundary faces       !
!   first, boundary region by boundary region, then inside faces, then         !
!   also accordin to indices of cells surrounding a face (c1 and c2).          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(grid_type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: s, c, c1, c2
  integer, allocatable :: new_c(:)
  integer, allocatable :: old_c(:)
  integer, allocatable :: i_work_1(:)
  integer, allocatable :: i_work_2(:,:)
  real,    allocatable :: criteria(:,:)
!==============================================================================!

  allocate(i_work_1(   grid % n_cells));     i_work_1(:)   = 0
  allocate(i_work_2(8, grid % n_cells));     i_work_2(:,:) = 0
  allocate(new_c   (   grid % n_cells));     new_c(:)      = 0
  allocate(old_c   (   grid % n_cells));     old_c(:)      = 0
  allocate(criteria(   grid % n_cells, 3));  criteria(:,:) = 0.0

  !------------------------!
  !   Store cells' nodes   !
  !------------------------!
  do c = 1, grid % n_cells
    i_work_1(c)      = grid % cells_n_nodes(c)
    i_work_2(1:8, c) = grid % cells_n(1:8, c)
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
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)
    grid % faces_c(1, s) = new_c(c1)
    if(c2 > 0) then
      grid % faces_c(2, s) = new_c(c2)
    end if
  end do

  !---------------------------!
  !   Do the actual sorting   !
  !---------------------------!
  do c = 1, grid % n_cells
    grid % cells_n_nodes(new_c(c)) = i_work_1(c)
    grid % cells_n(1:8, new_c(c))  = i_work_2(1:8, c)
  end do
  call Sort_Mod_Real_By_Index(grid % xc   (1), new_c(1), grid % n_cells)
  call Sort_Mod_Real_By_Index(grid % yc   (1), new_c(1), grid % n_cells)
  call Sort_Mod_Real_By_Index(grid % zc   (1), new_c(1), grid % n_cells)
  call Sort_Mod_Real_By_Index(grid % vol  (1), new_c(1), grid % n_cells)
  call Sort_Mod_Real_By_Index(grid % delta(1), new_c(1), grid % n_cells)

  end subroutine
