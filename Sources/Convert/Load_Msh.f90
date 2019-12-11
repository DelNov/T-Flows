!==============================================================================!
  subroutine Load_Msh(grid, mesh_city)
!------------------------------------------------------------------------------!
!   Reads the Gmsh file format.                                                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
  use Grid_Mod, only: Grid_Type,  &
                      Grid_Mod_Print_Bnd_Cond_List
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: mesh_city
!-----------------------------------[Locals]-----------------------------------!
  integer, parameter   :: MSH_TRI   = 2
  integer, parameter   :: MSH_QUAD  = 3
  integer, parameter   :: MSH_TETRA = 4
  integer, parameter   :: MSH_HEXA  = 5
  integer, parameter   :: MSH_PENTA = 6
  integer, parameter   :: MSH_PYRA  = 7
  character(len=130)   :: name_in
  integer              :: n_sect, n_elem, n_blocks, n_bnd_sect, n_grps, n_memb
  integer              :: i, j, c, dim, p_tag, s_tag, n_tags, type, fu
  integer              :: run, s_tag_max, n_e_0d, n_e_1d, n_e_2d, n_e_3d
  integer, allocatable :: n(:), new(:)
  integer, allocatable :: phys_tags(:), p_tag_corr(:), n_bnd_cells(:)
  character(len=80), allocatable :: phys_names(:)
!==============================================================================!

  call File_Mod_Set_Name(name_in, extension='.msh')
  call File_Mod_Open_File_For_Reading(name_in, fu)

  !----------------------------------------------!
  !                                              !
  !   Read number of nodes, cells, blocks, ...   !
  !                                              !
  !----------------------------------------------!

  !-------------------------------------------------!
  !   Read number of blocks and boundary sections   !
  !-------------------------------------------------!
  n_blocks   = 0
  n_bnd_sect = 0
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$PhysicalNames') exit
  end do
  call File_Mod_Read_Line(fu)
  read(line % tokens(1), *) n_sect
  allocate(phys_names(n_sect))
  allocate(p_tag_corr(n_sect))
  do i = 1, n_sect
    call File_Mod_Read_Line(fu)
    read(line % tokens(2), *) j
    if(line % tokens(1) .eq. '2') n_bnd_sect = n_bnd_sect + 1
    if(line % tokens(1) .eq. '3') n_blocks   = n_blocks   + 1
    read(line % tokens(2), *) j  ! section number; neglect
    if(line % tokens(1) .eq. '2') then
      read(line % tokens(3), *) phys_names(n_bnd_sect)
      p_tag_corr(j) = n_bnd_sect
    end if
  end do

  !--------------------------!
  !   Read number of nodes   !
  !--------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Nodes') exit
  end do
  call File_Mod_Read_Line(fu)
  read(line % tokens(4), *) grid % n_nodes  ! both 2 and 4 store number of nodes
  print *,'# Number of nodes: ', grid % n_nodes

  !--------------------------------------!
  !   Read number of elements (0D - 3D)  !
  !--------------------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Elements') exit
  end do
  call File_Mod_Read_Line(fu)
  read(line % tokens(4), *) n_elem  ! both 2 and 4 store number of elements
  allocate(new(n_elem))
  new(:) = 0

  !--------------------------------------!
  !   Read info on boundary conditions   !
  !--------------------------------------!
  do run = 1, 2  ! in the first run find max index
    if(run .eq. 1) s_tag_max = 0

    rewind(fu)
    do
      call File_Mod_Read_Line(fu)
      if(line % tokens(1) .eq. '$Entities') exit
    end do
    call File_Mod_Read_Line(fu)
    read(line % tokens(1), *) n_e_0d  ! number of 0D entities (points)
    read(line % tokens(2), *) n_e_1d  ! number of 1D entities (lines)
    read(line % tokens(3), *) n_e_2d  ! number of 2D entities (faces)
    read(line % tokens(4), *) n_e_3d  ! number of 3D entities (volumes)

    ! Skip 0D and 1D info
    do i = 1, n_e_0d + n_e_1d
      call File_Mod_Read_Line(fu)
    end do
    ! Analyze 2D data
    do i = 1, n_e_2d
      call File_Mod_Read_Line(fu)
      read(line % tokens(1), *) s_tag   ! surface tag
      read(line % tokens(8), *) n_tags  ! this should be one!  check some day
      read(line % tokens(9), *) p_tag   ! physcal tag
      if(n_tags .eq. 1) then
        if(run .eq. 1) s_tag_max = max(s_tag_max, s_tag)
        if(run .eq. 2) then
          phys_tags(s_tag) = p_tag_corr(p_tag)
        end if
      end if
      if(n_tags > 1) then
        print *, '# ERROR in Load_Msh @ s_tag: ', s_tag
        print *, '# More than one boundary condition per entity - not allowed!'
        stop
      end if
    end do
    if(run .eq. 1) then
      allocate(phys_tags(s_tag_max))
      phys_tags(:) = -1
    end if
  end do  ! next run

  !----------------------------------------!
  !   Count the inner and boundary cells   !
  !----------------------------------------!
  grid % n_bnd_cells = 0
  grid % n_cells     = 0
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Elements') exit
  end do
  call File_Mod_Read_Line(fu)
  read(line % tokens(1),*) n_grps
  do i = 1, n_grps
    call File_Mod_Read_Line(fu)
    read(line % tokens(1), *) dim     ! dimension of the element
    read(line % tokens(2), *) s_tag   ! element tag
    read(line % tokens(4), *) n_memb  ! number of members in the group
    do j = 1, n_memb
      call File_Mod_Read_Line(fu)
      read(line % tokens(1), *) c     ! Gmsh cell number
      if(dim .eq. 2) then
        grid % n_bnd_cells = grid % n_bnd_cells + 1
        new(c) = -grid % n_bnd_cells
      end if
      if(dim .eq. 3) then
        grid % n_cells = grid % n_cells + 1
        new(c) = grid % n_cells
      end if
    end do
  end do

  ! These five lines are coppied from Load_Neu
  print '(a38,i9)', '# Total number of nodes:             ', grid % n_nodes
  print '(a38,i9)', '# Total number of cells:             ', grid % n_cells
  print '(a38,i9)', '# Total number of blocks:            ', n_blocks
  print '(a38,i9)', '# Total number of boundary sections: ', n_bnd_sect
  print '(a38,i9)', '# Total number of boundary cells:    ', grid % n_bnd_cells

  !--------------------------------------------!
  !                                            !
  !   Allocate memory for Grid_Mod variables   !
  !                                            !
  !--------------------------------------------!
  call Allocate_Memory(grid)

  !---------------------------------------------------!
  !   Read boundary conditions for individual cells   !
  !---------------------------------------------------!
  allocate(n_bnd_cells(n_bnd_sect))
  n_bnd_cells(:) = 0

  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Elements') exit
  end do
  call File_Mod_Read_Line(fu)
  read(line % tokens(1),*) n_grps
  do i = 1, n_grps
    call File_Mod_Read_Line(fu)
    read(line % tokens(1), *) dim     ! dimension of the element
    read(line % tokens(2), *) s_tag   ! element tag
    read(line % tokens(4), *) n_memb  ! number of members in the group
    do j = 1, n_memb
      call File_Mod_Read_Line(fu)
      read(line % tokens(1), *) c     ! Gmsh cell number
      if(dim .eq. 2) then
        grid % bnd_cond % color( new(c) ) = phys_tags(s_tag)
        n_bnd_cells(phys_tags(s_tag)) = n_bnd_cells(phys_tags(s_tag)) + 1
      end if
    end do
  end do

  do i = 1, n_bnd_sect
    print '(a, i2, i6)', ' # Boundary sells in section: ', i, n_bnd_cells(i)
  end do

  !--------------------------------!
  !   Read the nodal coordinates   !
  !--------------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Nodes') exit
  end do
  call File_Mod_Read_Line(fu)
  read(line % tokens(1),*) n_grps    ! fetch number of groups
  do i = 1, n_grps
    call File_Mod_Read_Line(fu)
    read(line % tokens(4),*) n_memb  ! fetch number of members
    allocate(n(n_memb))
    do j = 1, n_memb                 ! fetch all node numbers
      call File_Mod_Read_Line(fu)
      read(line % tokens(1),*) n(j)
    end do
    do j = 1, n_memb
      call File_Mod_Read_Line(fu)    ! read node coordinates
      read(line % tokens(1),*) grid % xn(n(j))
      read(line % tokens(2),*) grid % yn(n(j))
      read(line % tokens(3),*) grid % zn(n(j))
    end do
    deallocate(n)
  end do

  !-----------------------------!
  !   Read nodes of each cell   !
  !-----------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Elements') exit
  end do
  call File_Mod_Read_Line(fu)
  read(line % tokens(1),*) n_grps
  do i = 1, n_grps
    call File_Mod_Read_Line(fu)
    read(line % tokens(1),*) dim     ! dimension of the element
    read(line % tokens(3),*) type    ! element type
    read(line % tokens(4),*) n_memb  ! number of members in the group
    do j = 1, n_memb
      call File_Mod_Read_Line(fu)

      ! Boundary cell, hopefully
      if(dim .eq. 2) then
        if(type .eq. MSH_TRI) then
          read(line % tokens(1), *) c       ! Gmsh cell number
          c = new(c)                        ! use T-Flows numbering
          grid % cells_n_nodes(c) = 3
          read(line % tokens(2), *) grid % cells_n(1, c)
          read(line % tokens(3), *) grid % cells_n(2, c)
          read(line % tokens(4), *) grid % cells_n(3, c)
        end if
        if(type .eq. MSH_QUAD) then
          read(line % tokens(1), *) c       ! Gmsh cell number
          c = new(c)                        ! use T-Flows numbering
          grid % cells_n_nodes(c) = 4
          read(line % tokens(2), *) grid % cells_n(1, c)
          read(line % tokens(3), *) grid % cells_n(2, c)
          read(line % tokens(4), *) grid % cells_n(3, c)
          read(line % tokens(5), *) grid % cells_n(4, c)
        end if
      end if

      ! Inside cells
      if(dim .eq. 3) then
        if(type .eq. MSH_TETRA) then
          read(line % tokens(1), *) c       ! Gmsh cell number
          c = new(c)                        ! use T-Flows numbering
          grid % cells_n_nodes(c) = 4
          read(line % tokens(2), *) grid % cells_n(1, c)
          read(line % tokens(3), *) grid % cells_n(2, c)
          read(line % tokens(4), *) grid % cells_n(3, c)
          read(line % tokens(5), *) grid % cells_n(4, c)
        end if
        if(type .eq. MSH_PENTA) then
          read(line % tokens(1), *) c       ! Gmsh cell number
          c = new(c)                        ! use T-Flows numbering
          grid % cells_n_nodes(c) = 6
          read(line % tokens(2), *) grid % cells_n(1, c)
          read(line % tokens(3), *) grid % cells_n(2, c)
          read(line % tokens(4), *) grid % cells_n(3, c)
          read(line % tokens(5), *) grid % cells_n(4, c)
          read(line % tokens(6), *) grid % cells_n(5, c)
          read(line % tokens(7), *) grid % cells_n(6, c)
        end if
        if(type .eq. MSH_HEXA) then
          read(line % tokens(1), *) c       ! Gmsh cell number
          c = new(c)                        ! use T-Flows numbering
          grid % cells_n_nodes(c) = 8
          read(line % tokens(2), *) grid % cells_n(1, c)
          read(line % tokens(3), *) grid % cells_n(2, c)
          read(line % tokens(4), *) grid % cells_n(4, c)
          read(line % tokens(5), *) grid % cells_n(3, c)
          read(line % tokens(6), *) grid % cells_n(5, c)
          read(line % tokens(7), *) grid % cells_n(6, c)
          read(line % tokens(8), *) grid % cells_n(8, c)
          read(line % tokens(9), *) grid % cells_n(7, c)
        end if
        if(type .eq. MSH_PYRA) then
          print *, '# ERROR: Pyramid cells not implemented yet!'
          stop
        end if
      end if
    end do
  end do

  !-----------------------!
  !                       !
  !   Set material name   !
  !                       !
  !-----------------------!
  grid % material % name = "AIR"

  !-------------------------!
  !   Boundary conditions   !
  !-------------------------!
  grid % n_bnd_cond = n_bnd_sect
  allocate(grid % bnd_cond % name(n_bnd_sect))

  do i = 1, n_bnd_sect
    grid % bnd_cond % name(i) = phys_names(i)
    call To_Upper_Case(grid % bnd_cond % name(i))
  end do

  !------------------------------------!
  !   Pring boundary conditions info   !
  !------------------------------------!
  call Grid_Mod_Print_Bnd_Cond_List(grid)

  close(fu)

  end subroutine
