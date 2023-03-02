!==============================================================================!
  subroutine Load_Gmsh(Convert, Grid, file_name)
!------------------------------------------------------------------------------!
!   Reads the Gmsh file format.                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert
  type(Grid_Type)     :: Grid
  character(SL)       :: file_name
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MSH_TRI   = 2
  integer, parameter :: MSH_QUAD  = 3
  integer, parameter :: MSH_TETRA = 4
  integer, parameter :: MSH_HEXA  = 5
  integer, parameter :: MSH_WEDGE = 6
  integer, parameter :: MSH_PYRA  = 7
!-----------------------------------[Locals]-----------------------------------!
  integer                    :: n_sect, n_elem, n_blocks, n_bnd_sect, n_grps
  integer                    :: n_memb, n_tags, n_crvs, n_nods, error
  integer                    :: i, j, k, c, dim, p_tag, s_tag, type, fu
  integer                    :: run, s_tag_max, n_e_0d, n_e_1d, n_e_2d, n_e_3d
  integer, allocatable       :: n(:), new(:)
  integer, allocatable       :: phys_tags(:), p_tag_corr(:), n_bnd_cells(:)
  character(SL), allocatable :: phys_names(:)
  logical                    :: ascii                 ! is file in ascii format?
  integer                    :: pos
  integer                    :: pos_meshformat     =  12
  integer                    :: pos_physicalnames  = -1
  integer                    :: pos_entities       = -1
  integer                    :: pos_nodes          = -1
  integer                    :: pos_elements       = -1
  integer(1)                 :: byte(0:3)
!==============================================================================!

  call Profiler % Start('Load_Gmsh')

  ! Open the file in binary mode, because it just might be
  call File % Open_For_Reading_Binary(file_name, fu)

  !-----------------------------------------------------------!
  !   A very rudimentary way to find beginnings of sections   !
  !-----------------------------------------------------------!
  pos = 0
  do
    read(fu, end=2) byte(0);  pos = pos + 1
    if(char(byte(0)) .eq. '$') then
      read(fu, end=2) byte(1);  pos = pos + 1
      read(fu, end=2) byte(2);  pos = pos + 1
      read(fu, end=2) byte(3);  pos = pos + 1
      if(char(byte(1)) .eq. 'E' .and.  &
         char(byte(2)) .eq. 'n' .and.  &
         char(byte(3)) .eq. 'd') then
        do i = 1, MAX_TOKENS*2
          read(fu, end=2) byte(0);  pos = pos + 1
          if(byte(0) .eq. 10) exit
        end do
        Line % whole = ''
        do i = 1, MAX_TOKENS*2
          read(fu, end=2) byte(0);  pos = pos + 1
          if(byte(0) .eq. 10) exit
          if(byte(0) .ne. 13) Line % whole(i:i) = char(byte(0))
        end do
        if(Line % whole .eq. '$PhysicalNames') pos_physicalnames = pos
        if(Line % whole .eq. '$Entities')      pos_entities      = pos
        if(Line % whole .eq. '$Nodes')         pos_nodes         = pos
        if(Line % whole .eq. '$Elements')      pos_elements      = pos
      end if
    end if
  end do
2 continue

  ! Error trap
  if(pos_physicalnames .eq. -1) then
    call Message % Error(60,                                              &
      "This is bad.  PhysicalNames don't seem to be defined in the "  //  &
      ".msh file.  Maybe you forgot to define boundary conditions  "  //  &
      "(called physical groups) in Gmsh?",                                &
      file=__FILE__, line=__LINE__)
  end if
  print *, '# Broswed the file in binary format and read ', pos, ' bytes'

  !----------------------------------------!
  !   Gmsh can't handle polyhedral grids   !
  !----------------------------------------!
  Grid % polyhedral = .false.

  !------------------------------!
  !   Check format fo the file   !
  !------------------------------!
  rewind(fu)
  call File % Read_Line(fu)
  Assert(Line % tokens(1) == '$MeshFormat')
  call File % Read_Line(fu)
  if(Line % tokens(1) .ne. '4.1') then
    call Message % Error(60,                                            &
             " Only Gmsh .msh files in version 4.1 are supported!"  //  &
             " \n I can't continue with this file, exiting!",           &
             file=__FILE__, line=__LINE__)
  end if

  ! My guess is that second tokens says it is binary (1) or not (0)
  ascii = .true.
  if(Line % tokens(2) .eq. '1') ascii = .false.

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
# ifdef __INTEL_COMPILER
  error = fseek(fu, pos_physicalnames, 0)
# else
  call fseek(fu, pos_physicalnames, 0)
# endif
  call File % Read_Line(fu)
  read(Line % tokens(1), *) n_sect
  allocate(phys_names(n_sect))
  allocate(p_tag_corr(n_sect * 65536))  ! allocate more than needed because
                                        ! it's all very messy in .msh files
  do i = 1, n_sect
    call File % Read_Line(fu)
    read(Line % tokens(2), *) j
    if(Line % tokens(1) .eq. '2') n_bnd_sect = n_bnd_sect + 1
    if(Line % tokens(1) .eq. '3') n_blocks   = n_blocks   + 1
    read(Line % tokens(2), *) j  ! section number; neglect
    if(Line % tokens(1) .eq. '2') then
      read(Line % tokens(3), *) phys_names(n_bnd_sect)
      p_tag_corr(j) = n_bnd_sect
    end if
  end do

  !--------------------------!
  !                          !
  !   Read number of nodes   !
  !                          !
  !--------------------------!
  rewind(fu)
# ifdef __INTEL_COMPILER
  error = fseek(fu, pos_nodes, 0)
# else
  call fseek(fu, pos_nodes, 0)
# endif
  if(ascii) then
    call File % Read_Line(fu)
    read(Line % tokens(4), *) Grid % n_nodes  ! 2 and 4 store number of nodes
  else
    call File % Read_Binary_Int8_Array(fu, 4)
    Grid % n_nodes = int8_array(4)
  end if
  print *,'# Number of nodes: ', Grid % n_nodes

  !--------------------------------------!
  !                                      !
  !   Read number of elements (0D - 3D)  !
  !                                      !
  !--------------------------------------!
  rewind(fu)
# ifdef __INTEL_COMPILER
  error = fseek(fu, pos_elements, 0)
# else
  call fseek(fu, pos_elements, 0)
# endif
  if(ascii) then
    call File % Read_Line(fu)
    read(Line % tokens(4), *) n_elem  ! both 2 and 4 store number of elements
  else
    call File % Read_Binary_Int8_Array(fu, 4)
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
#   ifdef __INTEL_COMPILER
    error = fseek(fu, pos_entities, 0)
#   else
    call fseek(fu, pos_entities, 0)
#   endif
    if(ascii) then
      call File % Read_Line(fu)
      read(Line % tokens(1), *) n_e_0d  ! number of 0D entities (points)
      read(Line % tokens(2), *) n_e_1d  ! number of 1D entities (lines)
      read(Line % tokens(3), *) n_e_2d  ! number of 2D entities (faces)
      read(Line % tokens(4), *) n_e_3d  ! number of 3D entities (volumes)
    else
      call File % Read_Binary_Int8_Array(fu, 4)
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
        call File % Read_Line(fu)
      end do
    else
      do i = 1, n_e_0d
        ! Node's tag
        call File % Read_Binary_Int4_Array (fu, 1)
        ! Node's coordinates
        call File % Read_Binary_Real8_Array(fu, 3)
        ! Number of physical tags (it is assumed to be zero, to check maybe?)
        call File % Read_Binary_Int8_Array (fu, 1)
      end do
    end if

    !--------------------------!
    !   Skip 1D (curve) info   !
    !--------------------------!
    if(ascii) then
      do i = 1, n_e_1d
        call File % Read_Line(fu)
      end do
    else
      do i = 1, n_e_1d
        ! Curve's tag
        call File % Read_Binary_Int4_Array (fu, 1)
        ! Bounding box coordinates
        call File % Read_Binary_Real8_Array(fu, 6)
        ! Number of physical tags (it is assumed to be zero, to check maybe?)
        call File % Read_Binary_Int8_Array (fu, 1)
        ! Number of bounding points (assumed to be two, a check one day?)
        call File % Read_Binary_Int8_Array (fu, 1)
        ! Points one and two
        call File % Read_Binary_Int4_Array (fu, 2)
      end do
    end if

    !-------------------------------!
    !   Analyze 2D (surface) data   !
    !-------------------------------!
    do i = 1, n_e_2d

      if(ascii) then
        call File % Read_Line(fu)
        read(Line % tokens(1), *) s_tag          ! surface tag
        read(Line % tokens(8), *) n_tags         ! for me this was 1 or 0
        do j = 1, n_tags                         ! browse physical tags ...
          read(Line % tokens(8+j), *) p_tag      ! ... and read them
        end do
        read(Line % tokens(9+n_tags), *) n_crvs  ! number of bounding curves
      else
        ! Surface's tag
        call File % Read_Binary_Int4_Array (fu, 1)
        s_tag = int4_array(1)
        ! Bounding box coordinates
        call File % Read_Binary_Real8_Array(fu, 6)
        ! Number of physical tags
        call File % Read_Binary_Int8_Array (fu, 1)
        n_tags = int8_array(1)
        do j = 1, n_tags  ! read the physical tags you have
          call File % Read_Binary_Int4_Array (fu, 1)
          p_tag = int4_array(1)
        end do
        ! Number of bounding curves
        call File % Read_Binary_Int8_Array (fu, 1)
        n_crvs = int8_array(1)
        ! Read the bounding curves
        call File % Read_Binary_Int4_Array (fu, n_crvs)
      end if
      if(n_tags .eq. 1) then
        if(run .eq. 1) s_tag_max = max(s_tag_max, s_tag)
        if(run .eq. 2) then
          phys_tags(s_tag) = p_tag_corr(p_tag)
        end if
      end if
      if(n_tags > 1) then
        call Message % Error(50,                                        &
               " More than one boundary condition per entity. \n "  //  &
               " It is not supported in this verion of T-Flows!",       &
               file=__FILE__, line=__LINE__)
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
  Grid % n_bnd_cells = 0
  Grid % n_cells     = 0
  rewind(fu)
# ifdef __INTEL_COMPILER
  error = fseek(fu, pos_elements, 0)
# else
  call fseek(fu, pos_elements, 0)
# endif

  !-------------------------------!
  !   Read the number of groups   !
  !-------------------------------!
  if(ascii) then
    call File % Read_Line(fu)
    read(Line % tokens(1),*) n_grps
  else
    call File % Read_Binary_Int8_Array(fu, 4)
    n_grps = int8_array(1)
  end if

  !-------------------------------------------------------------!
  !   Browse through groups to count boundary and inner cells   !
  !-------------------------------------------------------------!
  do i = 1, n_grps

    ! Read dim, s_tag, type and n_memb <--= this is what you actually need!
    if(ascii) then
      call File % Read_Line(fu)
      read(Line % tokens(1), *) dim     ! dimension of the element
      read(Line % tokens(2), *) s_tag   ! element tag
      read(Line % tokens(3), *) type    ! element type
      read(Line % tokens(4), *) n_memb  ! number of members in the group
    else
      call File % Read_Binary_Int4_Array(fu, 3)
      dim   = int4_array(1)  ! dimension of the element
      s_tag = int4_array(2)  ! element tag
      type  = int4_array(3)  ! element type
      call File % Read_Binary_Int8_Array(fu, 1)
      n_memb = int8_array(1)  ! number of members in the group
    end if

    ! Read cell number and cell's nodes <--= this is just to carry on
    do j = 1, n_memb
      if(ascii) then
        call File % Read_Line(fu)
        read(Line % tokens(1), *) c     ! Gmsh cell number
      else
        ! Element tag
        call File % Read_Binary_Int8_Array(fu, 1)
        c = int8_array(1)
        ! Node tags
        if(type .eq. MSH_TRI)   call File % Read_Binary_Int8_Array(fu, 3)
        if(type .eq. MSH_QUAD)  call File % Read_Binary_Int8_Array(fu, 4)
        if(type .eq. MSH_TETRA) call File % Read_Binary_Int8_Array(fu, 4)
        if(type .eq. MSH_HEXA)  call File % Read_Binary_Int8_Array(fu, 8)
        if(type .eq. MSH_WEDGE) call File % Read_Binary_Int8_Array(fu, 6)
        if(type .eq. MSH_PYRA)  call File % Read_Binary_Int8_Array(fu, 5)
      end if
      if(dim .eq. 2) then
        Grid % n_bnd_cells = Grid % n_bnd_cells + 1
        new(c) = -Grid % n_bnd_cells
      end if
      if(dim .eq. 3) then
        Grid % n_cells = Grid % n_cells + 1
        new(c) = Grid % n_cells
      end if
    end do
  end do    ! n_grps

  ! These five lines are coppied from Load_Neu
  print '(a38,i9)', '# Total number of nodes:             ', Grid % n_nodes
  print '(a38,i9)', '# Total number of cells:             ', Grid % n_cells
  print '(a38,i9)', '# Total number of blocks:            ', n_blocks
  print '(a38,i9)', '# Total number of boundary sections: ', n_bnd_sect
  print '(a38,i9)', '# Total number of boundary cells:    ', Grid % n_bnd_cells

  !--------------------------------------------!
  !                                            !
  !   Allocate memory for Grid_Mod variables   !
  !                                            !
  !--------------------------------------------!
  call Convert % Allocate_Memory(Grid)

  !---------------------------------------------------!
  !                                                   !
  !   Read boundary conditions for individual cells   !
  !                                                   !
  !---------------------------------------------------!
  allocate(n_bnd_cells(n_bnd_sect))
  n_bnd_cells(:) = 0

  rewind(fu)
# ifdef __INTEL_COMPILER
  error = fseek(fu, pos_elements, 0)
# else
  call fseek(fu, pos_elements, 0)
# endif

  !-------------------------------!
  !   Read the number of groups   !
  !-------------------------------!
  if(ascii) then
    call File % Read_Line(fu)
    read(Line % tokens(1),*) n_grps
  else
    call File % Read_Binary_Int8_Array(fu, 4)
    n_grps = int8_array(1)
  end if

  !----------------------------------------------------------------!
  !   Browse through groups, read and extract more detailed info   !
  !----------------------------------------------------------------!
  do i = 1, n_grps

    ! Read dim, s_tag, type and n_memb
    if(ascii) then
      call File % Read_Line(fu)
      read(Line % tokens(1), *) dim     ! dimension of the element
      read(Line % tokens(2), *) s_tag   ! element tag
      read(Line % tokens(3), *) type    ! element type
      read(Line % tokens(4), *) n_memb  ! number of members in the group
    else
      call File % Read_Binary_Int4_Array(fu, 3)
      dim   = int4_array(1)  ! dimension of the element
      s_tag = int4_array(2)  ! element tag
      type  = int4_array(3)  ! element type
      call File % Read_Binary_Int8_Array(fu, 1)
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
        call File % Read_Line(fu)
        read(Line % tokens(1), *) c  ! fetch Gmsh cell number
        c = new(c)                   ! use T-Flows numbering

        Grid % cells_n_nodes(c) = n_nods
        call Adjust_First_Dim(n_nods, Grid % cells_n)
        do k = 1, n_nods
          read(Line % tokens(k+1), *) Grid % cells_n(k, c)
        end do

      else  ! it is in binary format

        ! Element tag
        call File % Read_Binary_Int8_Array(fu, n_nods+1)
        c = int8_array(1)  ! fetch Gmsh cell number
        c = new(c)         ! use T-Flows numbering

        Grid % cells_n_nodes(c) = n_nods
        call Adjust_First_Dim(n_nods, Grid % cells_n)
        do k = 1, n_nods
          Grid % cells_n(k, c) = int8_array(k+1)
        end do

      end if

      if(dim .eq. 2) then
        Grid % region % at_cell(c) = phys_tags(s_tag)
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
# ifdef __INTEL_COMPILER
  error = fseek(fu, pos_nodes, 0)
# else
  call fseek(fu, pos_nodes, 0)
# endif

  !-------------------------------!
  !   Read the number of groups   !
  !-------------------------------!
  if(ascii) then
    call File % Read_Line(fu)
    read(Line % tokens(1),*) n_grps
  else
    call File % Read_Binary_Int8_Array(fu, 4)
    n_grps = int8_array(1)
  end if

  !-------------------------------------------------------!
  !   Browse through groups and fetch nodal coordinates   !
  !-------------------------------------------------------!
  do i = 1, n_grps
    if(ascii) then
      call File % Read_Line(fu)
      read(Line % tokens(4),*) n_memb  ! fetch number of members
    else
      call File % Read_Binary_Int4_Array(fu, 3)
      call File % Read_Binary_Int8_Array(fu, 1)
      n_memb = int8_array(1)
    end if
    allocate(n(n_memb))

    ! Fetch all node numbers in the group
    if(ascii) then
      do j = 1, n_memb
        call File % Read_Line(fu)
        read(Line % tokens(1),*) n(j)
      end do
    else
      do j = 1, n_memb                 ! fetch all node numbers
        call File % Read_Binary_Int8_Array(fu, 1)
        n(j) = int8_array(1)
      end do
    end if

    ! Fetch all node coordinates in the group
    if(ascii) then
      do j = 1, n_memb
        call File % Read_Line(fu)    ! read node coordinates
        read(Line % tokens(1),*) Grid % xn(n(j))
        read(Line % tokens(2),*) Grid % yn(n(j))
        read(Line % tokens(3),*) Grid % zn(n(j))
      end do
    else
      do j = 1, n_memb
        call File % Read_Binary_Real8_Array(fu, 3)
        Grid % xn(n(j)) = real8_array(1)
        Grid % yn(n(j)) = real8_array(2)
        Grid % zn(n(j)) = real8_array(3)
      end do
    end if
    deallocate(n)
  end do

  !----------------------------------!
  !                                  !
  !   Copy boundary condition info   !
  !                                  !
  !----------------------------------!
  call Grid % Allocate_Regions(n_bnd_sect)

  do i = 1, n_bnd_sect
    Grid % region % name(i) = phys_names(i)
    call String % To_Upper_Case(Grid % region % name(i))
  end do

  !------------------------------------!
  !                                    !
  !   Print boundary conditions info   !
  !                                    !
  !------------------------------------!
  call Grid % Print_Regions_List()

  close(fu)

  call Profiler % Stop('Load_Gmsh')

  end subroutine
