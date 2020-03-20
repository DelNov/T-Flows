!==============================================================================!
  subroutine Grid_Mod_Load_Cns(grid, this_proc)
!------------------------------------------------------------------------------!
!   Reads: .cns file.                                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: this_proc  ! needed if called from Processor
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n, s, lev, fu
  character(len=80) :: name_in
!==============================================================================!

  !-------------------------------!
  !                               !
  !     Read the file with the    !
  !   connections between cells   !
  !                               !
  !-------------------------------!
  call File_Mod_Set_Name(name_in, processor=this_proc, extension='.cns')
  call File_Mod_Open_File_For_Reading_Binary(name_in, fu, this_proc)

  !-----------------------------------------------!
  !   Number of cells, boundary cells and faces   !
  !-----------------------------------------------!
  read(fu) grid % n_nodes
  read(fu) grid % n_cells              ! number of cells including buffer
  read(fu) grid % n_bnd_cells          ! number of boundary cells
  read(fu) grid % n_faces              ! number of faces (with buffer faces)
  read(fu) grid % comm % n_buff_cells  ! number of buffer faces/cells
  read(fu) grid % n_bnd_cond           ! number of boundary conditions
  read(fu) grid % n_levels             ! number of multigrid levels

  ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells, grid % n_bnd_cells)
  call Grid_Mod_Allocate_Faces(grid, grid % n_faces)

  ! Boundary conditions' keys
  allocate(grid % bnd_cond % name(grid % n_bnd_cond))
  allocate(grid % bnd_cond % type(grid % n_bnd_cond))

  !-------------------!
  !   Material name   !
  !-------------------!
  read(fu) grid % material % name

  !------------------------------!
  !   Boundary conditions list   !
  !------------------------------!
  do n = 1, grid % n_bnd_cond
    read(fu) grid % bnd_cond % name(n)
  end do

  !-----------!
  !   Cells   !  (including buffer cells)
  !-----------!

  ! Number of nodes for each cell
  read(fu) (grid % cells_n_nodes(c), c = 1, grid % n_cells)

  ! Cells' nodes
  read(fu) ((grid % cells_n(n,c),              &
             n = 1, grid % cells_n_nodes(c)),  &
             c = 1, grid % n_cells)

  ! Cells' processor ids
  read(fu) (grid % comm % cell_proc(c), c =  1,  grid % n_cells)
  read(fu) (grid % comm % cell_proc(c), c = -1, -grid % n_bnd_cells, -1)

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  read(fu) (grid % faces_n_nodes(s), s = 1, grid % n_faces)

  ! Faces' nodes
  read(fu) ((grid % faces_n(n,s),              &
             n = 1, grid % faces_n_nodes(s)),  &
             s = 1, grid % n_faces)

  ! Faces' cells
  read(fu) ((grid % faces_c(c,s), c = 1, 2), s = 1, grid % n_faces)

  !--------------!
  !   Boundary   !
  !--------------!

  ! Physical boundary cells
  allocate (grid % bnd_cond % color(-grid % n_bnd_cells-1:-1))
  read(fu) (grid % bnd_cond % color(c), c = -1,-grid % n_bnd_cells, -1)

  call Grid_Mod_Bnd_Cond_Ranges(grid)

  !----------------------!
  !   Multigrid levels   !
  !----------------------!
  do lev = 1, grid % n_levels
    read(fu)  grid % level(lev) % n_cells
    read(fu)  grid % level(lev) % n_faces
  end do
  call Grid_Mod_Allocate_Levels(grid)
  do lev = 1, grid % n_levels
    read(fu) (grid % level(lev) % cell(c),      c=1,grid % n_cells)
    read(fu) (grid % level(lev) % face(s),      s=1,grid % n_faces)
    read(fu) (grid % level(lev) % coarser_c(c), c=1,grid % level(lev) % n_cells)
    read(fu) (grid % level(lev) % faces_c(1,s), s=1,grid % level(lev) % n_faces)
    read(fu) (grid % level(lev) % faces_c(2,s), s=1,grid % level(lev) % n_faces)
  end do
  call Grid_Mod_Check_Levels(grid)

  close(fu)

  end subroutine
