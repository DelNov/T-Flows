!==============================================================================!
  subroutine Grid_Mod_Load_Cns(grid, this_proc, domain)
!------------------------------------------------------------------------------!
!   Reads: .cns file.                                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  integer           :: this_proc  ! needed if called from Processor
  integer, optional :: domain
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, s, lev, fu
  character(SL) :: name_in
!==============================================================================!

  !-------------------------------!
  !                               !
  !     Read the file with the    !
  !   connections between cells   !
  !                               !
  !-------------------------------!
  call File_Mod_Set_Name(name_in, processor=this_proc, extension='.cns',  &
                         domain=domain)
  call File_Mod_Open_File_For_Reading_Binary(name_in, fu, this_proc)

  !-----------------------------------------------!
  !   Number of cells, boundary cells and faces   !
  !-----------------------------------------------!
  read(fu) grid % n_nodes
  read(fu) grid % n_cells              ! number of cells including buffer
  read(fu) grid % n_bnd_cells          ! number of boundary cells
  read(fu) grid % n_faces              ! number of faces (with buffer faces)
  read(fu) grid % n_shadows            ! number of shadow faces
  read(fu) grid % n_bnd_cond           ! number of boundary conditions

  ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells, grid % n_bnd_cells)
  call Grid_Mod_Allocate_Faces(grid, grid % n_faces, grid % n_shadows)

  ! Boundary conditions' keys
  ! (Go from zero for faces which are not at the boundary)
  allocate(grid % bnd_cond % name(0 : grid % n_bnd_cond + 3))
  allocate(grid % bnd_cond % type(0 : grid % n_bnd_cond + 3))

  !-----------------!
  !   Domain name   !
  !-----------------!
  read(fu) grid % name

  !------------------------------!
  !   Boundary conditions list   !
  !------------------------------!
  do n = 1, grid % n_bnd_cond
    read(fu) grid % bnd_cond % name(n)
  end do
  ! The last three are reserved for perodicity
  ! (Assigned in Grid_Mod_Find_Periodic_Faces)
  grid % bnd_cond % name(grid % n_bnd_cond + 1) = 'PERIODIC_X'
  grid % bnd_cond % name(grid % n_bnd_cond + 2) = 'PERIODIC_Y'
  grid % bnd_cond % name(grid % n_bnd_cond + 3) = 'PERIODIC_Z'

  !--------------------------!
  !   Nodes global numbers   !
  !--------------------------!
  read(fu) (grid % comm % node_glo(n), n = 1, grid % n_nodes)

  !-----------!
  !   Cells   !  (including buffer cells)
  !-----------!

  ! Number of nodes for each cell
  read(fu) (grid % cells_n_nodes(c), c = -grid % n_bnd_cells, grid % n_cells)

  ! Cells' nodes
  read(fu) ((grid % cells_n(n,c),              &
             n = 1, grid % cells_n_nodes(c)),  &
             c = -grid % n_bnd_cells, grid % n_cells)

  ! Cells' processor ids
  read(fu) (grid % comm % cell_proc(c), c = -grid % n_bnd_cells, grid % n_cells)

  ! Cells' global indices
  read(fu) (grid % comm % cell_glo(c), c = -grid % n_bnd_cells, grid % n_cells)

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  read(fu) (grid % faces_n_nodes(s), s = 1, grid % n_faces + grid % n_shadows)

  ! Faces' nodes
  read(fu) ((grid % faces_n(n,s),              &
             n = 1, grid % faces_n_nodes(s)),  &
             s = 1, grid % n_faces + grid % n_shadows)

  ! Faces' cells
  read(fu) ((grid % faces_c(c,s), c = 1, 2), s = 1, grid % n_faces  &
                                                  + grid % n_shadows)

  ! Faces' shadows
  read(fu) (grid % faces_s(s), s = 1, grid % n_faces + grid % n_shadows)

  ! Faces' global numbers
  read(fu) (grid % comm % face_glo(s), s = 1, grid % n_faces + grid % n_shadows)

  !--------------!
  !   Boundary   !
  !--------------!

  ! Physical boundary cells (and all the faces)
  ! (This opens the oportunity to store bounary condition info in ...
  !  ... the faces thus ridding us of the "if(c2 < 0) then" checks)
  allocate (grid % bnd_cond % color(-grid % n_bnd_cells-1:grid % n_faces))
  read(fu) (grid % bnd_cond % color(c), c = -grid % n_bnd_cells, -1)

  call Grid_Mod_Bnd_Cond_Ranges(grid)

  close(fu)

  end subroutine
