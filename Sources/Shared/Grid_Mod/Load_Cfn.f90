!==============================================================================!
  subroutine Grid_Mod_Load_Cfn(grid, this_proc, domain)
!------------------------------------------------------------------------------!
!   Reads: .cfn file.                                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  integer           :: this_proc  ! needed if called from Processor
  integer, optional :: domain
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, c1, c2, s, n, ss, sr, fu
  character(SL) :: name_in
!==============================================================================!

  !-------------------------------!
  !                               !
  !     Read the file with the    !
  !   connections between cells   !
  !                               !
  !-------------------------------!
  call File_Mod_Set_Name(name_in,              &
                         processor=this_proc,  &
                         extension='.cfn',     &
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

  !-------------------------------------!
  !   Does grid have polyhedral cells   !
  !-------------------------------------!
  read(fu) grid % polyhedral

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
  ! and used for inlet copy boundary condition.  Don't delete these thinking
  ! they are useless.  They are assigned in Grid_Mod_Calculate_Face_Geometry
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

  ! Error trap for number of nodes for each cell
  do c = -grid % n_bnd_cells, grid % n_cells
    if(c .ne. 0) then
      if(grid % cells_n_nodes(c) .eq. 0) then
        print *, '# ERROR: Number of nodes is zero at cell:', c
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end if
  end do

  ! Cells' nodes
  read(fu) ((grid % cells_n(n, c),                     &
             n = 1, abs(grid % cells_n_nodes(c))),     &
             c = -grid % n_bnd_cells, grid % n_cells)

  ! Error trap for cells' nodes
  do c = -grid % n_bnd_cells, grid % n_cells
    do n = 1, abs(grid % cells_n_nodes(c))
      if(grid % cells_n(n, c) .eq. 0) then
        print *, '# ERROR: Node index is zero at cell:', c
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end do
  end do

  ! Number of faces for each cell
  read(fu) (grid % cells_n_faces(c), c = -grid % n_bnd_cells, grid % n_cells)

  ! Error trap for number of faces for each cell
  do c = -grid % n_bnd_cells, grid % n_cells
    if(c .ne. 0) then
      if(grid % cells_n_faces(c) .eq. 0) then
        print *, '# ERROR: Number of faces is zero at cell:', c
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end if
  end do

  ! Cells' faces
  read(fu) ((grid % cells_f(s, c),             &
             s = 1, grid % cells_n_faces(c)),  &
             c = -grid % n_bnd_cells, grid % n_cells)

  ! Error trap for cells' faces
  do c = -grid % n_bnd_cells, grid % n_cells
    do s = 1, grid % cells_n_faces(c)
      if(grid % cells_f(s, c) .eq. 0) then
        print *, '# ERROR: Face index is zero at cell:', c
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end do
  end do

  ! Cells' processor ids
  read(fu) (grid % comm % cell_proc(c), c = -grid % n_bnd_cells, grid % n_cells)

  ! Cells' global indices
  read(fu) (grid % comm % cell_glo(c), c = -grid % n_bnd_cells, grid % n_cells)

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  read(fu) (grid % faces_n_nodes(s), s = 1, grid % n_faces + grid % n_shadows)

  ! Error trap for number of nodes for each face
  do s = 1, grid % n_faces + grid % n_shadows
    if(grid % faces_n_nodes(s) .eq. 0) then
      print *, '# ERROR: Number of nodes is zero at face:', s
      print *, '# This error is critical.  Exiting!'
      call Comm_Mod_End
      stop
    end if
  end do

  ! Faces' nodes
  read(fu) ((grid % faces_n(n, s),             &
             n = 1, grid % faces_n_nodes(s)),  &
             s = 1, grid % n_faces + grid % n_shadows)

  ! Error trap for faces' nodes
  do s = 1, grid % n_faces + grid % n_shadows
    do n = 1, grid % faces_n_nodes(s)
      if(grid % faces_n(n, s) .eq. 0) then
        print *, '# ERROR: Node index is zero at face:', s
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end do
  end do

  ! Faces' cells
  read(fu) ((grid % faces_c(c, s), c = 1, 2), s = 1, grid % n_faces  &
                                                   + grid % n_shadows)

  ! Error trap for faces' cells
  do s = 1, grid % n_faces + grid % n_shadows
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    ! Check only if least one cell is in this processor
    ! (Meaning it is not a face entirelly in the buffer)
    if(grid % comm % cell_proc(c1) .eq. this_proc .or.  &
       grid % comm % cell_proc(c2) .eq. this_proc) then
      if( .not. (c1.eq.0 .and. c2.eq.0) ) then
        if(grid % faces_c(1, s) .eq. 0) then
          print *, '# ERROR: Cell one is zero at face:', s, c1, c2
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
        if(grid % faces_c(2, s) .eq. 0) then
          print *, '# ERROR: Cell two is zero at face:', s
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end if
    end if
  end do

  ! Faces' shadows
  read(fu) (grid % faces_s(s), s = 1, grid % n_faces + grid % n_shadows)

  ! Error trap for shadows
  do ss = grid % n_faces + 1, grid % n_faces + grid % n_shadows
    sr = grid % faces_s(ss)  ! real face from shadow data
    if(sr .eq. 0) then
      print *, '# ERROR: Shadow faces points to zero face'
      print *, '# This error is critical.  Exiting!'
      call Comm_Mod_End
      stop
    end if
    if(grid % faces_s(sr) .ne. ss) then
      print *, '# ERROR: Real and shadow faces do not point to each other'
      print *, '# This error is critical.  Exiting!'
      call Comm_Mod_End
      stop
    end if
  end do

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
