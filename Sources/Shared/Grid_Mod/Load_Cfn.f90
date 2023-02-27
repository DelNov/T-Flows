!==============================================================================!
  subroutine Load_Cfn(Grid, this_proc, domain)
!------------------------------------------------------------------------------!
!   Reads: .cfn file.                                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  integer, intent(in) :: this_proc  ! needed if called from Processor
  integer, optional   :: domain
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, c1, c2, s, n, ss, sr, fu, real_prec
  character(SL) :: name_in
!==============================================================================!

  !-------------------------------!
  !     Read the file with the    !
  !   connections between cells   !
  !-------------------------------!
  call File % Set_Name(name_in,              &
                       processor=this_proc,  &
                       extension='.cfn',     &
                       domain=domain)
  call File % Open_For_Reading_Binary(name_in, fu, this_proc)

  !-------------------------!
  !   Read real precision   !
  !-------------------------!
  read(fu) real_prec

  if(real_prec .ne. RP) then
    if(RP .eq. 8 .and. real_prec .eq. 4) then
      call Message % Error(64,                                                 &
                       'Input files were saved in single precision, but '   // &
                       'this program is compiled in double. Decide in   '   // &
                       'which precision you want to work, recompile the '   // &
                       'programs (option REAL=double/single in make), '     // &
                       'convert or generate the meshes again and re-run.',     &
                       one_proc = .true.)
    else
      call Message % Error(64,                                                 &
                       'Input files were saved in double precision, but '   // &
                       'this program is compiled in single. Decide in   '   // &
                       'which precision you want to work, recompile the '   // &
                       'programs (option REAL=double/single in make), '     // &
                       'convert or generate the meshes again and re-run.',     &
                       one_proc = .true.)
    end if
  end if

  !-----------------------------------------------!
  !   Number of cells, boundary cells and faces   !
  !-----------------------------------------------!
  read(fu) Grid % n_nodes
  read(fu) Grid % n_cells              ! number of cells including buffer
  read(fu) Grid % n_bnd_cells          ! number of boundary cells
  read(fu) Grid % n_faces              ! number of faces (with buffer faces)
  read(fu) Grid % n_shadows            ! number of shadow faces
  read(fu) Grid % n_bnd_regions        ! number of boundary conditions

  !-------------------------------------!
  !   Does grid have polyhedral cells   !
  !-------------------------------------!
  read(fu) Grid % polyhedral

  ! Allocate memory =--> carefull, there is no checking!
  call Grid % Allocate_Nodes  (Grid % n_nodes)
  call Grid % Allocate_Cells  (Grid % n_cells, Grid % n_bnd_cells)
  call Grid % Allocate_Faces  (Grid % n_faces, Grid % n_shadows)
  call Grid % Allocate_Regions(Grid % n_bnd_regions)

  !-----------------!
  !   Domain name   !
  !-----------------!
  read(fu) Grid % name

  !------------------------------!
  !   Boundary conditions list   !
  !------------------------------!
  do n = Boundary_Regions()
    read(fu) Grid % region % name(n)
  end do

  !--------------------------!
  !   Nodes global numbers   !
  !--------------------------!
  read(fu) (Grid % Comm % node_glo(n), n = 1, Grid % n_nodes)

  !-----------!
  !   Cells   !  (including buffer cells)
  !-----------!

  ! Number of nodes for each cell
  read(fu) (Grid % cells_n_nodes(c), c = -Grid % n_bnd_cells, Grid % n_cells)

  ! Error trap for number of nodes for each cell
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(c .ne. 0) then
      if(Grid % cells_n_nodes(c) .eq. 0) then
        print *, '# ERROR: Number of nodes is zero at cell:', c
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end if
  end do

  ! Cells' nodes
  do c = -Grid % n_bnd_cells, Grid % n_cells
    call Adjust_First_Dim(abs(Grid % cells_n_nodes(c)),  &
                          Grid % cells_n)
    do n = 1, abs(Grid % cells_n_nodes(c))
      read(fu) Grid % cells_n(n, c)
    end do
  end do

  ! Error trap for cells' nodes
  do c = -Grid % n_bnd_cells, Grid % n_cells
    do n = 1, abs(Grid % cells_n_nodes(c))
      if(Grid % cells_n(n, c) .eq. 0) then
        print *, '# ERROR: Node index is zero at cell:', c
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end do
  end do

  ! Number of faces for each cell
  read(fu) (Grid % cells_n_faces(c), c = -Grid % n_bnd_cells, Grid % n_cells)

  ! Error trap for number of faces for each cell
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(c .ne. 0) then
      if(Grid % cells_n_faces(c) .eq. 0) then
        print *, '# ERROR: Number of faces is zero at cell:', c
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end if
  end do

  ! Cells' faces
  do c = -Grid % n_bnd_cells, Grid % n_cells
    call Adjust_First_Dim(Grid % cells_n_faces(c), Grid % cells_f)
    do s = 1, Grid % cells_n_faces(c)
      read(fu) Grid % cells_f(s, c)
    end do
  end do

  ! Error trap for cells' faces
  do c = -Grid % n_bnd_cells, Grid % n_cells
    do s = 1, Grid % cells_n_faces(c)
      if(Grid % cells_f(s, c) .eq. 0) then
        print *, '# ERROR: Face index is zero at cell:', c
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end do
  end do

  ! Cells' processor ids
  read(fu) (Grid % Comm % cell_proc(c), c = -Grid % n_bnd_cells, Grid % n_cells)

  ! Cells' global indices
  read(fu) (Grid % Comm % cell_glo(c), c = -Grid % n_bnd_cells, Grid % n_cells)

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  read(fu) (Grid % faces_n_nodes(s), s = 1, Grid % n_faces + Grid % n_shadows)

  ! Error trap for number of nodes for each face
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % faces_n_nodes(s) .eq. 0) then
      print *, '# ERROR: Number of nodes is zero at face:', s
      print *, '# This error is critical.  Exiting!'
      call Comm_Mod_End
      stop
    end if
  end do

  ! Faces' nodes
  do s = 1, Grid % n_faces + Grid % n_shadows
    call Adjust_First_Dim(Grid % faces_n_nodes(s), Grid % faces_n)
    do n = 1, Grid % faces_n_nodes(s)
      read(fu) Grid % faces_n(n, s)
    end do
  end do

  ! Error trap for faces' nodes
  do s = 1, Grid % n_faces + Grid % n_shadows
    do n = 1, Grid % faces_n_nodes(s)
      if(Grid % faces_n(n, s) .eq. 0) then
        print *, '# ERROR: Node index is zero at face:', s
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end do
  end do

  ! Faces' cells
  read(fu) ((Grid % faces_c(c, s), c = 1, 2), s = 1, Grid % n_faces  &
                                                   + Grid % n_shadows)

  ! Error trap for faces' cells
  do s = 1, Grid % n_faces + Grid % n_shadows
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    ! Check only if least one cell is in this processor
    ! (Meaning it is not a face entirelly in the buffer)
    if(Grid % Comm % cell_proc(c1) .eq. this_proc .or.  &
       Grid % Comm % cell_proc(c2) .eq. this_proc) then
      if( .not. (c1.eq.0 .and. c2.eq.0) ) then
        if(Grid % faces_c(1, s) .eq. 0) then
          print *, '# ERROR: Cell one is zero at face:', s, c1, c2
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
        if(Grid % faces_c(2, s) .eq. 0) then
          print *, '# ERROR: Cell two is zero at face:', s
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end if
    end if
  end do

  ! Faces' shadows
  read(fu) (Grid % faces_s(s), s = 1, Grid % n_faces + Grid % n_shadows)

  ! Error trap for shadows
  do ss = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows
    sr = Grid % faces_s(ss)  ! real face from shadow data
    if(sr .eq. 0) then
      print *, '# ERROR: Shadow faces points to zero face'
      print *, '# This error is critical.  Exiting!'
      call Comm_Mod_End
      stop
    end if
    if(Grid % faces_s(sr) .ne. ss) then
      print *, '# ERROR: Real and shadow faces do not point to each other'
      print *, '# This error is critical.  Exiting!'
      call Comm_Mod_End
      stop
    end if
  end do

  !--------------!
  !   Boundary   !
  !--------------!

  ! Physical boundary cells (and all the faces)
  ! (This opens the oportunity to store bounary condition info in ...
  !  ... the faces thus ridding us of the "if(c2 < 0) then" checks)
  allocate (Grid % region % at_cell(-Grid % n_bnd_cells-1:Grid % n_faces))
  read(fu) (Grid % region % at_cell(c), c = -Grid % n_bnd_cells, -1)

  call Grid % Determine_Regions_Ranges()

  close(fu)

  end subroutine
