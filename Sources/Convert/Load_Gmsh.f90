!==============================================================================!
  subroutine Load_Gmsh(grid, file_name)
!------------------------------------------------------------------------------!
!   Reads the Gmsh file format.                                                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  character(SL)   :: file_name
!-----------------------------------[Locals]-----------------------------------!
  integer                    :: n_sect, n_elem, n_blocks, n_bnd_sect, n_grps
  integer                    :: n_memb, n_tags, n_crvs, n_nods
  integer                    :: i, j, c, dim, p_tag, s_tag, type, fu
  integer                    :: run, s_tag_max, n_e_0d, n_e_1d, n_e_2d, n_e_3d
  integer, allocatable       :: n(:), new(:)
  integer, allocatable       :: phys_tags(:), p_tag_corr(:), n_bnd_cells(:)
  character(SL), allocatable :: phys_names(:)
  logical                    :: ascii                 ! is file in ascii format?
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MSH_TRI   = 2
  integer, parameter :: MSH_QUAD  = 3
  integer, parameter :: MSH_TETRA = 4
  integer, parameter :: MSH_HEXA  = 5
  integer, parameter :: MSH_WEDGE = 6
  integer, parameter :: MSH_PYRA  = 7
!==============================================================================!

  ! Open the file in binary mode, because it just might be
  call File_Mod_Open_File_For_Reading_Binary(file_name, fu)

  !----------------------------------------!
  !   Gmsh can't handle polyhedral grids   !
  !----------------------------------------!
  grid % polyhedral = .false.

  !------------------------------!
  !   Check format fo the file   !
  !------------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$MeshFormat') exit
  end do
  call File_Mod_Read_Line(fu)
  if(line % tokens(1) .ne. '4.1') then
    print *, '# ERROR in Load_Gmsh: files in version 4.1 are supported!'
    print *, '# This error is criticial.  Exiting!'
    stop
  end if

  ! My guess is that second tokens says it is binary (1) or not (0)
  ascii = .true.
  if(line % tokens(2) .eq. '1') ascii = .false.
  ! Line which follows contains some crap, but who cares?

  !----------------------------------------------!
  !                                              !
  !                                              !
  !   Read number of nodes, cells, blocks, ...   !
  !                                              !
  !                                              !
  !----------------------------------------------!

  !-------------------------------------------------!
  !                                                 !
  !   Read number of blocks and boundary sections   !
  !                                                 !
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
  allocate(p_tag_corr(n_sect * 128))  ! allocate more than needed because
                                      ! it's all very messy in .msh files
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
  !                          !
  !   Read number of nodes   !
  !                          !
  !--------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Nodes') exit
  end do
  if(ascii) then
    call File_Mod_Read_Line(fu)
    read(line % tokens(4), *) grid % n_nodes  ! 2 and 4 store number of nodes
  else
    call File_Mod_Read_Binary_Int8_Array(fu, 4)
    grid % n_nodes = int8_array(4)
  end if
  print *,'# Number of nodes: ', grid % n_nodes

  !--------------------------------------!
  !                                      !
  !   Read number of elements (0D - 3D)  !
  !                                      !
  !--------------------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Elements') exit
  end do
  if(ascii) then
    call File_Mod_Read_Line(fu)
    read(line % tokens(4), *) n_elem  ! both 2 and 4 store number of elements
  else
    call File_Mod_Read_Binary_Int8_Array(fu, 4)
    n_elem = int8_array(4)
  end if
  allocate(new(n_elem))
  new(:) = 0

  !--------------------------------------!
  !                                      !
  !   Read info on boundary conditions   !
  !                                      !
  !--------------------------------------!
  do run = 1, 2  ! in the first run find max index
    if(run .eq. 1) s_tag_max = 0

    rewind(fu)
    do
      call File_Mod_Read_Line(fu)
      if(line % tokens(1) .eq. '$Entities') exit
    end do
    if(ascii) then
      call File_Mod_Read_Line(fu)
      read(line % tokens(1), *) n_e_0d  ! number of 0D entities (points)
      read(line % tokens(2), *) n_e_1d  ! number of 1D entities (lines)
      read(line % tokens(3), *) n_e_2d  ! number of 2D entities (faces)
      read(line % tokens(4), *) n_e_3d  ! number of 3D entities (volumes)
    else
      call File_Mod_Read_Binary_Int8_Array(fu, 4)
      n_e_0d = int8_array(1)  ! number of 0D entities (points)
      n_e_1d = int8_array(2)  ! number of 1D entities (lines)
      n_e_2d = int8_array(3)  ! number of 2D entities (faces)
      n_e_3d = int8_array(4)  ! number of 3D entities (volumes)
    end if

    !--------------------------!
    !   Skip 0D (point) info   !
    !--------------------------!
    if(ascii) then
      do i = 1, n_e_0d
        call File_Mod_Read_Line(fu)
      end do
    else
      do i = 1, n_e_0d
        ! Node's tag
        call File_Mod_Read_Binary_Int4_Array (fu, 1)
        ! Node's coordinates
        call File_Mod_Read_Binary_Real8_Array(fu, 3)
        ! Number of physical tags (it is assumed to be zero, to check maybe?)
        call File_Mod_Read_Binary_Int8_Array (fu, 1)
      end do
    end if

    !--------------------------!
    !   Skip 1D (curve) info   !
    !--------------------------!
    if(ascii) then
      do i = 1, n_e_1d
        call File_Mod_Read_Line(fu)
      end do
    else
      do i = 1, n_e_1d
        ! Curve's tag
        call File_Mod_Read_Binary_Int4_Array (fu, 1)
        ! Bounding box coordinates
        call File_Mod_Read_Binary_Real8_Array(fu, 6)
        ! Number of physical tags (it is assumed to be zero, to check maybe?)
        call File_Mod_Read_Binary_Int8_Array (fu, 1)
        ! Number of bounding points (assumed to be two, a check one day?)
        call File_Mod_Read_Binary_Int8_Array (fu, 1)
        ! Points one and two
        call File_Mod_Read_Binary_Int4_Array (fu, 2)
      end do
    end if

    !-------------------------------!
    !   Analyze 2D (surface) data   !
    !-------------------------------!
    do i = 1, n_e_2d

      if(ascii) then
        call File_Mod_Read_Line(fu)
        read(line % tokens(1), *) s_tag          ! surface tag
        read(line % tokens(8), *) n_tags         ! for me this was 1 or 0
        do j = 1, n_tags                         ! browse physical tags ...
          read(line % tokens(8+j), *) p_tag      ! ... and read them
        end do
        read(line % tokens(9+n_tags), *) n_crvs  ! number of bounding curves
      else
        ! Surface's tag
        call File_Mod_Read_Binary_Int4_Array (fu, 1)
        s_tag = int4_array(1)
        ! Bounding box coordinates
        call File_Mod_Read_Binary_Real8_Array(fu, 6)
        ! Number of physical tags
        call File_Mod_Read_Binary_Int8_Array (fu, 1)
        n_tags = int8_array(1)
        do j = 1, n_tags  ! read the physical tags you have
          call File_Mod_Read_Binary_Int4_Array (fu, 1)
          p_tag = int4_array(1)
        end do
        ! Number of bounding curves
        call File_Mod_Read_Binary_Int8_Array (fu, 1)
        n_crvs = int8_array(1)
        ! Read the bounding curves
        call File_Mod_Read_Binary_Int4_Array (fu, n_crvs)
      end if
      if(n_tags .eq. 1) then
        if(run .eq. 1) s_tag_max = max(s_tag_max, s_tag)
        if(run .eq. 2) then
          phys_tags(s_tag) = p_tag_corr(p_tag)
        end if
      end if
      if(n_tags > 1) then
        print *, '# ERROR in Load_Gmsh @ s_tag: ', s_tag
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
  !                                        !
  !   Count the inner and boundary cells   !
  !                                        !
  !----------------------------------------!
  grid % n_bnd_cells = 0
  grid % n_cells     = 0
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Elements') exit
  end do

  !-------------------------------!
  !   Read the number of groups   !
  !-------------------------------!
  if(ascii) then
    call File_Mod_Read_Line(fu)
    read(line % tokens(1),*) n_grps
  else
    call File_Mod_Read_Binary_Int8_Array(fu, 4)
    n_grps = int8_array(1)
  end if

  !-------------------------------------------------------------!
  !   Browse through groups to count boundary and inner cells   !
  !-------------------------------------------------------------!
  do i = 1, n_grps

    ! Read dim, s_tag, type and n_memb <--= this is what you actually need!
    if(ascii) then
      call File_Mod_Read_Line(fu)
      read(line % tokens(1), *) dim     ! dimension of the element
      read(line % tokens(2), *) s_tag   ! element tag
      read(line % tokens(3), *) type    ! element type
      read(line % tokens(4), *) n_memb  ! number of members in the group
    else
      call File_Mod_Read_Binary_Int4_Array(fu, 3)
      dim   = int4_array(1)  ! dimension of the element
      s_tag = int4_array(2)  ! element tag
      type  = int4_array(3)  ! element type
      call File_Mod_Read_Binary_Int8_Array(fu, 1)
      n_memb = int8_array(1)  ! number of members in the group
    end if

    ! Read cell number and cell's nodes <--= this is just to carry on
    do j = 1, n_memb
      if(ascii) then
        call File_Mod_Read_Line(fu)
        read(line % tokens(1), *) c     ! Gmsh cell number
      else
        ! Element tag
        call File_Mod_Read_Binary_Int8_Array(fu, 1)
        c = int8_array(1)
        ! Node tags
        if(type .eq. MSH_TRI)   call File_Mod_Read_Binary_Int8_Array(fu, 3)
        if(type .eq. MSH_QUAD)  call File_Mod_Read_Binary_Int8_Array(fu, 4)
        if(type .eq. MSH_TETRA) call File_Mod_Read_Binary_Int8_Array(fu, 4)
        if(type .eq. MSH_HEXA)  call File_Mod_Read_Binary_Int8_Array(fu, 8)
        if(type .eq. MSH_WEDGE) call File_Mod_Read_Binary_Int8_Array(fu, 6)
        if(type .eq. MSH_PYRA)  call File_Mod_Read_Binary_Int8_Array(fu, 5)
      end if
      if(dim .eq. 2) then
        grid % n_bnd_cells = grid % n_bnd_cells + 1
        new(c) = -grid % n_bnd_cells
      end if
      if(dim .eq. 3) then
        grid % n_cells = grid % n_cells + 1
        new(c) = grid % n_cells
      end if
    end do
  end do    ! n_grps

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
  !                                                   !
  !   Read boundary conditions for individual cells   !
  !                                                   !
  !---------------------------------------------------!
  allocate(n_bnd_cells(n_bnd_sect))
  n_bnd_cells(:) = 0

  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Elements') exit
  end do

  !-------------------------------!
  !   Read the number of groups   !
  !-------------------------------!
  if(ascii) then
    call File_Mod_Read_Line(fu)
    read(line % tokens(1),*) n_grps
  else
    call File_Mod_Read_Binary_Int8_Array(fu, 4)
    n_grps = int8_array(1)
  end if

  !----------------------------------------------------------------!
  !   Browse through groups, read and extract more detailed info   !
  !----------------------------------------------------------------!
  do i = 1, n_grps

    ! Read dim, s_tag, type and n_memb
    if(ascii) then
      call File_Mod_Read_Line(fu)
      read(line % tokens(1), *) dim     ! dimension of the element
      read(line % tokens(2), *) s_tag   ! element tag
      read(line % tokens(3), *) type    ! element type
      read(line % tokens(4), *) n_memb  ! number of members in the group
    else
      call File_Mod_Read_Binary_Int4_Array(fu, 3)
      dim   = int4_array(1)  ! dimension of the element
      s_tag = int4_array(2)  ! element tag
      type  = int4_array(3)  ! element type
      call File_Mod_Read_Binary_Int8_Array(fu, 1)
      n_memb = int8_array(1)  ! number of members in the group
    end if

    ! Treat different cell types now
    if(type .eq. MSH_TRI)   n_nods = 3
    if(type .eq. MSH_QUAD)  n_nods = 4
    if(type .eq. MSH_TETRA) n_nods = 4
    if(type .eq. MSH_WEDGE) n_nods = 6
    if(type .eq. MSH_HEXA)  n_nods = 8
    if(type .eq. MSH_PYRA)  n_nods = 5

    ! Read cell number and cell's nodes
    do j = 1, n_memb
      if(ascii) then
        call File_Mod_Read_Line(fu)
        read(line % tokens(1), *) c  ! fetch Gmsh cell number
        c = new(c)                   ! use T-Flows numbering

        grid % cells_n_nodes(c) = n_nods
        read(line % tokens(2:n_nods+1), *) grid % cells_n(1:n_nods, c)

      else  ! it is in binary format

        ! Element tag
        call File_Mod_Read_Binary_Int8_Array(fu, n_nods+1)
        c = int8_array(1)  ! fetch Gmsh cell number
        c = new(c)         ! use T-Flows numbering

        grid % cells_n_nodes(c) = n_nods
        grid % cells_n(1:n_nods, c) = int8_array(2:n_nods+1)

      end if

      if(dim .eq. 2) then
        grid % bnd_cond % color(c) = phys_tags(s_tag)
        n_bnd_cells(phys_tags(s_tag)) = n_bnd_cells(phys_tags(s_tag)) + 1
      end if
    end do

  end do  ! n_grps

  do i = 1, n_bnd_sect
    print '(a, i2, i7)', ' # Boundary cells in section: ', i, n_bnd_cells(i)
  end do

  !--------------------------------!
  !                                !
  !   Read the nodal coordinates   !
  !                                !
  !--------------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. '$Nodes') exit
  end do

  !-------------------------------!
  !   Read the number of groups   !
  !-------------------------------!
  if(ascii) then
    call File_Mod_Read_Line(fu)
    read(line % tokens(1),*) n_grps
  else
    call File_Mod_Read_Binary_Int8_Array(fu, 4)
    n_grps = int8_array(1)
  end if

  !-------------------------------------------------------!
  !   Browse through groups and fetch nodal coordinates   !
  !-------------------------------------------------------!
  do i = 1, n_grps
    if(ascii) then
      call File_Mod_Read_Line(fu)
      read(line % tokens(4),*) n_memb  ! fetch number of members
    else
      call File_Mod_Read_Binary_Int4_Array(fu, 3)
      call File_Mod_Read_Binary_Int8_Array(fu, 1)
      n_memb = int8_array(1)
    end if
    allocate(n(n_memb))

    ! Fetch all node numbers in the group
    if(ascii) then
      do j = 1, n_memb
        call File_Mod_Read_Line(fu)
        read(line % tokens(1),*) n(j)
      end do
    else
      do j = 1, n_memb                 ! fetch all node numbers
        call File_Mod_Read_Binary_Int8_Array(fu, 1)
        n(j) = int8_array(1)
      end do
    end if

    ! Fetch all node coordinates in the group
    if(ascii) then
      do j = 1, n_memb
        call File_Mod_Read_Line(fu)    ! read node coordinates
        read(line % tokens(1),*) grid % xn(n(j))
        read(line % tokens(2),*) grid % yn(n(j))
        read(line % tokens(3),*) grid % zn(n(j))
      end do
    else
      do j = 1, n_memb
        call File_Mod_Read_Binary_Real8_Array(fu, 3)
        grid % xn(n(j)) = real8_array(1)
        grid % yn(n(j)) = real8_array(2)
        grid % zn(n(j)) = real8_array(3)
      end do
    end if
    deallocate(n)
  end do

  !----------------------------------!
  !                                  !
  !   Copy boundary condition info   !
  !                                  !
  !----------------------------------!
  grid % n_bnd_cond = n_bnd_sect
  allocate(grid % bnd_cond % name(n_bnd_sect))

  do i = 1, n_bnd_sect
    grid % bnd_cond % name(i) = phys_names(i)
    call To_Upper_Case(grid % bnd_cond % name(i))
  end do

  !------------------------------------!
  !                                    !
  !   Print boundary conditions info   !
  !                                    !
  !------------------------------------!
  call Grid_Mod_Print_Bnd_Cond_List(grid)

  close(fu)

  end subroutine
