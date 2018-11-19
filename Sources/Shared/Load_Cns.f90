!==============================================================================!
  subroutine Load_Cns(grid, this_proc)
!------------------------------------------------------------------------------!
!   Reads: .cns file.                                                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: this_proc  ! needed if called from Processor
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n, s
  character(len=80) :: name_in
!==============================================================================!

  !-------------------------------!
  !     Read the file with the    !
  !   connections between cells   !
  !-------------------------------!
  call Name_File(this_proc, name_in, '.cns')  

  open(9, file=name_in,form='unformatted')
  if(this_proc < 2) print *, '# Reading the file: ', name_in

  ! Number of cells, boundary cells and sides
  read(9) grid % n_nodes
  read(9) grid % n_cells              ! number of cells including buffer
  read(9) grid % n_bnd_cells          ! number of boundary cells
  read(9) grid % n_faces              ! number of faces (with buffer faces)
  read(9) grid % comm % n_buff_cells  ! number of buffer faces/cells
  read(9) grid % n_bnd_cond           ! number of boundary conditions

  ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells, grid % n_bnd_cells) 
  call Grid_Mod_Allocate_Faces(grid, grid % n_faces) 

  ! Boundary conditions' keys
  allocate(grid % bnd_cond % name(grid % n_bnd_cond))
  allocate(grid % bnd_cond % type(grid % n_bnd_cond))
  read(9) grid % material % name
  do n = 1, grid % n_bnd_cond
    read(9) grid % bnd_cond % name(n)
  end do

  ! Cells
  read(9) (grid % cells_n_nodes(c), c = 1, grid % n_cells)
  read(9) ((grid % cells_n(n,c),              &
            n = 1, grid % cells_n_nodes(c)),  &
            c = 1, grid % n_cells)

  ! Cells' processor ids
  read(9) (grid % comm % proces(c), c =  1, grid % n_cells)
  read(9) (grid % comm % proces(c), c = -1,-grid % n_bnd_cells,-1)

  ! Faces
  read(9) (grid % faces_n_nodes(s), s = 1, grid % n_faces)
  read(9) ((grid % faces_n(n,s),              &
            n = 1, grid % faces_n_nodes(s)),  &
            s = 1, grid % n_faces)
  read(9) (grid % faces_c(1,s), s = 1, grid % n_faces)
  read(9) (grid % faces_c(2,s), s = 1, grid % n_faces)

  ! Boundary cells
  allocate (grid % bnd_cond % color(-grid % n_bnd_cells-1:-1))
  read(9) (grid % bnd_cond % color(c), c = -1,-grid % n_bnd_cells, -1)

  call Grid_Mod_Bnd_Cond_Ranges(grid)

  ! Boundary copy cells
  allocate (grid % bnd_cond % copy_c(-grid % n_bnd_cells:-1))
  read(9) (grid % bnd_cond % copy_c(c), c = -1,-grid % n_bnd_cells, -1)

  read(9) grid % n_copy
  allocate (grid % bnd_cond % copy_s(2,grid % n_copy))
  read(9) (grid % bnd_cond % copy_s(1,s), s = 1,grid % n_copy)
  read(9) (grid % bnd_cond % copy_s(2,s), s = 1,grid % n_copy)

  close(9)

  end subroutine
