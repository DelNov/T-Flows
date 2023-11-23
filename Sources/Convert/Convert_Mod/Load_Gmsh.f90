!==============================================================================!
  subroutine Load_Gmsh(Convert, Grid, file_name)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to read grid files in the Gmsh file format
!>  and populate the provided Grid object with the necessary data.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initial setup: The subroutine starts by opening the specified Gmsh file  !
!     in binary mode for initial analysis.                                     !
!   * Section identification: It identifies important sections in the Gmsh     !
!     file, such as PhysicalNames, Entities, Nodes, and Elements, using a      !
!     buffered reading approach. (This is crucial for locating and parsing     !
!     different parts of the file efficiently.)                                !
!   * File format check: The subroutine checks the version of the Gmsh file    !
!     and whether it's in ASCII or binary format. This affects how the file    !
!     will be read and processed.                                              !
!   * Reading grid information:                                                !
!     - It reads the number of blocks and boundary sections in the grid.       !
!     - The subroutine counts the number of internal and boundary cells, a key !
!       step for allocating memory for the grid data.                          !
!     - It processes boundary conditions and entities, extracting data such    !
!       as surface dimensions, element types, and the number of members in     !
!       each group.                                                            !
!   * Memory allocation for grid: Allocates memory for the Grid object based   !
!     on the extracted information, preparing it for storing grid data.        !
!   * Reading cell and node data:                                              !
!     - It reads and stores cell information, including cell types (triangles, !
!       quadrilaterals, tetrahedra, etc.) and node indices.                    !
!     - The subroutine also reads and stores nodal coordinates.                !
!   * Boundary condition processing: Reads and processes boundary condition    !
!     information, translating Gmsh's physical group tags to T-Flows' format.  !
!   * Finalization: The subroutine closes the file, stops profiling, and       !
!     prepares the Grid object with all the necessary grid data for further    !
!     processing.                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert    !! parent class
  type(Grid_Type)     :: Grid       !! grid being converted
  character(SL)       :: file_name  !! file name
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MSH_TRI   = 2
  integer, parameter :: MSH_QUAD  = 3
  integer, parameter :: MSH_TETRA = 4
  integer, parameter :: MSH_HEXA  = 5
  integer, parameter :: MSH_WEDGE = 6
  integer, parameter :: MSH_PYRA  = 7
  integer, parameter :: L         = 1048576  ! buffer size
  integer, parameter :: W         =       8  ! buffer's window size
  integer, parameter :: P         =       5  ! pattern's length
!-----------------------------------[Locals]-----------------------------------!
  type(Pattern_Type)         :: PhysicalNames
  type(Pattern_Type)         :: Entities
  type(Pattern_Type)         :: Nodes
  type(Pattern_Type)         :: Elements
  integer                    :: pos_physicalnames  = -1
  integer                    :: pos_entities       = -1
  integer                    :: pos_nodes          = -1
  integer                    :: pos_elements       = -1
  integer                    :: n_sect, n_elem, n_blocks, n_bnd_sect, n_grps
  integer                    :: n_memb, n_tags, n_crvs, n_nods, error, ios
  integer                    :: i, j, k, c, dim, p_tag, s_tag, type, fu, tot
  integer                    :: run, s_tag_max, n_e_0d, n_e_1d, n_e_2d, n_e_3d
  integer                    :: f, e  ! buffer indices: first (f) and end (e)
  integer,       allocatable :: n(:), new(:)
  integer,       allocatable :: phys_tags(:), p_tag_corr(:), n_bnd_cells(:)
  character(SL), allocatable :: phys_names(:)
  logical                    :: ascii                 ! is file in ascii format?
  integer(1)                 :: buffer(L)
!------------------------------------------------------------------------------!
!   Buffered reading in search of patterns:                                    !
!                                                                              !
!   The way to handle this is to maintain a sliding window that keeps track    !
!   of the last few bytes of the previous buffer. When you load new buffer,    !
!   you combine the last few bytes of the old buffer with the new buffer.      !
!   Thus, even if pattern is split between buffers, you can still detect it.   !
!                                                                              !
!   Representation of a buffer and explanation of some local variables:        !
!                                                                              !
!   + is the leading window from the previous buffer                           !
!   o is the current buffer                                                    !
!                                                                              !
!   buffer:    ++++++++oooooooooooooooooooooooooooooooooooooooooooooooo        !
!   global:    f       s                                              e        !
!   local:     1       W                                              L        !
!                                                                              !
!   There is no point in searching for a pattern from s till e, since you      !
!   know that you won't be able to find it if it split by buffers.  Hence      !
!   you search from s - P (where P is the pattern lengt) till e - P.  The      !
!   trailing part of the buffer will be coppied to the next sweep              !
!==============================================================================!

  call Profiler % Start('Load_Gmsh')

  !-------------------------------------------------------!
  !                                                       !
  !   Open the file in binary mode for initial analyzis   !
  !                                                       !
  !-------------------------------------------------------!
  call File % Open_For_Reading_Binary(file_name, fu)

  call Profiler % Start('Load_Gmsh (find sections)')

  !-------------------!
  !                   !
  !   Find sections   !
  !                   !
  !-------------------!

  ! A chance that another five character pattern will occur
  ! in a binary file is 1/256^5 or 9.094947e-13 or 1.e-12)
  call PhysicalNames % Create_Pattern("$Phys")
  call Entities      % Create_Pattern("$Enti")
  call Nodes         % Create_Pattern("$Node")
  call Elements      % Create_Pattern("$Elem")

  !---------------------------------------------!
  !   Read the input file in windowed buffers   !
  !---------------------------------------------!
  tot = 0  ! number of items already found
  do

    ! Get the first position, or position before reading
    inquire(unit=fu, pos=f)  ! get the global position within the file
                             ! when you read the file for the first time,
                             ! it will be one, pos points to last+!

    ! Read data into the buffer, from the end of the window (W) till the end
    read(fu, iostat=ios) buffer(W+1:L)

    ! Get the position after reading
    inquire(unit=fu, pos=e)  ! get the global position within the file

    Assert(ios .le. 0)       ! check if an error occurred

    ! At this position, I have at my disposal buffer from local position
    ! W+1 till L to analyze, which corresponds to global positions s to e
    ! However, since I want to analyze from s - P till e - P, I should
    ! browse from W + 1 - P to L - P
    do i = W+1-P, L-P
      if( PhysicalNames % Match_Pattern(buffer(i)) ) then
        pos_physicalnames  = i+f-1-W + len("$PhysicalNames")
        tot = tot + 1
      end if
      if( Entities % Match_Pattern(buffer(i)) ) then
        pos_entities       = i+f-1-W + len("$Entities")
        tot = tot + 1
      end if
      if( Nodes % Match_Pattern(buffer(i)) ) then
        pos_nodes          = i+f-1-W + len("$Nodes")
        tot = tot + 1
      end if
      if( Elements % Match_Pattern(buffer(i)) ) then
        pos_elements       = i+f-1-W + len("$Elements")
        tot = tot + 1
      end if
      if(tot .eq. 4) goto 1
    end do

    ! Copy the last window back for the next buffer
    buffer(1:W) = buffer(L-W+1:L)

    if(ios < 0) exit  ! if end of file reached, exit
  end do

1 continue

  call Profiler % Stop('Load_Gmsh (find sections)')

  ! Error trap
  if(pos_physicalnames .eq. -1) then
    call Message % Error(60,                                              &
      "This is bad.  PhysicalNames don't seem to be defined in the "  //  &
      ".msh file.  Maybe you forgot to define boundary conditions  "  //  &
      "(called physical groups) in Gmsh?",                                &
      file=__FILE__, line=__LINE__)
  end if

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

  !------------------------------------------!
  !                                          !
  !   Once you know the format of the file   !
  !     close it and open again if ASCII     !
  !                                          !
  !------------------------------------------!
  if(ascii) then
    close(fu)
    call File % Open_For_Reading_Ascii(file_name, fu)
  end if

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
  error = 0
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
    read(fu, *) (Line % tokens(k), k = 1, 4)
    read(Line % tokens(4), *) Grid % n_nodes  ! 2 and 4 store number of nodes
  else
    call File % Read_Binary_Int8_Array(fu, 4)
    Grid % n_nodes = int(int8_array(4))
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
    read(fu, *) (Line % tokens(k), k = 1, 4)
    read(Line % tokens(4), *) n_elem  ! both 2 and 4 store number of elements
  else
    call File % Read_Binary_Int8_Array(fu, 4)
    n_elem = int(int8_array(4))
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
      read(fu, *) (Line % tokens(k), k = 1, 4)
      read(Line % tokens(1), *) n_e_0d  ! number of 0D entities (points)
      read(Line % tokens(2), *) n_e_1d  ! number of 1D entities (lines)
      read(Line % tokens(3), *) n_e_2d  ! number of 2D entities (faces)
      read(Line % tokens(4), *) n_e_3d  ! number of 3D entities (volumes)
    else
      call File % Read_Binary_Int8_Array(fu, 4)
      n_e_0d = int(int8_array(1))  ! number of 0D entities (points)
      n_e_1d = int(int8_array(2))  ! number of 1D entities (lines)
      n_e_2d = int(int8_array(3))  ! number of 2D entities (faces)
      n_e_3d = int(int8_array(4))  ! number of 3D entities (volumes)
    end if

    !--------------------------!
    !   Skip 0D (point) info   !
    !--------------------------!
    if(ascii) then
      do i = 1, n_e_0d
        read(fu, *) Line % whole
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
        read(fu, *) Line % whole
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
        n_tags = int(int8_array(1))
        do j = 1, n_tags  ! read the physical tags you have
          call File % Read_Binary_Int4_Array (fu, 1)
          p_tag = int4_array(1)
        end do
        ! Number of bounding curves
        call File % Read_Binary_Int8_Array (fu, 1)
        n_crvs = int(int8_array(1))
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
    read(fu, *) n_grps
  else
    call File % Read_Binary_Int8_Array(fu, 4)
    n_grps = int(int8_array(1))
  end if

  !-------------------------------------------------------------!
  !   Browse through groups to count boundary and inner cells   !
  !-------------------------------------------------------------!
  do i = 1, n_grps

    ! Read dim, s_tag, type and n_memb <--= this is what you actually need!
    if(ascii) then
      read(fu, *) (Line % tokens(k), k = 1, 4)
      read(Line % tokens(1), *) dim     ! dimension of the element
      read(Line % tokens(2), *) s_tag   ! element tag
      read(Line % tokens(3), *) type    ! element type
      read(Line % tokens(4), *) n_memb  ! number of members in the group
    else
      call File % Read_Binary_Int4_Array(fu, 3)
      dim   = int4_array(1)        ! dimension of the element
      s_tag = int4_array(2)        ! element tag
      type  = int4_array(3)        ! element type
      call File % Read_Binary_Int8_Array(fu, 1)
      n_memb = int(int8_array(1))  ! number of members in the group
    end if

    ! Read cell number and cell's nodes <--= this is just to carry on
    do j = 1, n_memb
      if(ascii) then
        read(fu, *) c
      else
        ! Element tag
        call File % Read_Binary_Int8_Array(fu, 1)
        c = int(int8_array(1))
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

  if(n_bnd_sect .eq. 0) then
    call Message % Error(60, 'No boundary sections (physical groups in  ' //  &
                             'Gmsh) have been defined.  Convert can''t '  //  &
                             'work with grids like that. '                //  &
                             '\n \n '                                     //  &
                             'Define physical groups in Gmsh and retry.',     &
                             file=__FILE__, line=__LINE__)
  end if

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
    read(fu, *) n_grps
  else
    call File % Read_Binary_Int8_Array(fu, 4)
    n_grps = int(int8_array(1))
  end if

  !----------------------------------------------------------------!
  !   Browse through groups, read and extract more detailed info   !
  !----------------------------------------------------------------!
  print '(a)', ' # Browse through groups, for more detailed info'
  do i = 1, n_grps  ! there are hundreds of these groups :-(

    ! Read dim, s_tag, type and n_memb
    if(ascii) then
      read(fu, *) (Line % tokens(k), k = 1, 4)
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
      n_memb = int(int8_array(1))  ! number of members in the group
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
        read(fu, *) (Line % tokens(k), k = 1, 1 + n_nods)
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
        c = int(int8_array(1))  ! fetch Gmsh cell number
        c = new(c)              ! use T-Flows numbering

        Grid % cells_n_nodes(c) = n_nods
        call Adjust_First_Dim(n_nods, Grid % cells_n)
        do k = 1, n_nods
          Grid % cells_n(k, c) = int(int8_array(k+1))
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
    read(fu, *) n_grps
  else
    call File % Read_Binary_Int8_Array(fu, 4)
    n_grps = int(int8_array(1))
  end if

  !-------------------------------------------------------!
  !   Browse through groups and fetch nodal coordinates   !
  !-------------------------------------------------------!
  do i = 1, n_grps
    if(ascii) then
      read(fu, *) (Line % tokens(k), k = 1, 4)
      read(Line % tokens(4),*) n_memb  ! fetch number of members
    else
      call File % Read_Binary_Int4_Array(fu, 3)
      call File % Read_Binary_Int8_Array(fu, 1)
      n_memb = int(int8_array(1))
    end if
    allocate(n(n_memb))

    ! Fetch all node numbers in the group
    if(ascii) then
      do j = 1, n_memb
        read(fu, *) n(j)
      end do
    else
      do j = 1, n_memb                 ! fetch all node numbers
        call File % Read_Binary_Int8_Array(fu, 1)
        n(j) = int(int8_array(1))
      end do
    end if

    ! Fetch all node coordinates in the group
    if(ascii) then
      do j = 1, n_memb
        read(fu, *) Grid % xn(n(j)), Grid % yn(n(j)), Grid % zn(n(j))
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
