!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Connections(base, block, sect, grid)
!------------------------------------------------------------------------------!
!   Read elements connection info for current sect
!
!   In case of PW, 2d section connection node has the same name as b.c. 
!   in ZoneBC -> not reliable
!   In case of Salome it is not true
!   Same with interfaces
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer         :: base, block, sect
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base_id       ! base index number
  integer              :: block_id      ! block index number
  integer              :: sect_id       ! element section index
  character(len=80)    :: sect_name     ! name of the Elements_t node
  integer              :: cell_type     ! types of elements in the section
  integer              :: first_cell    ! index of first element
  integer              :: last_cell     ! index of last element
  integer              :: parent_flag
  integer              :: error
  integer              :: n_nodes, loc, c, n, cnt, pos
  integer, allocatable :: cell_n(:,:)
  integer, allocatable :: mixed_cell_n(:), &
                          cell_n_hex(:,:), &
                          cell_n_pyr(:,:), &
                          cell_n_tet(:,:), &
                          cell_n_wed(:,:)
  integer              :: loc_hex, loc_pyr, loc_wed, loc_tet
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  sect_id  = sect

  ! Introduce some abbreviations
  sect_name   = cgns_base(base) % block(block) % section(sect) % name
  cell_type   = cgns_base(base) % block(block) % section(sect) % cell_type
  first_cell  = cgns_base(base) % block(block) % section(sect) % first_cell
  last_cell   = cgns_base(base) % block(block) % section(sect) % last_cell
  parent_flag = cgns_base(base) % block(block) % section(sect) % parent_flag

  ! Number of cells in this section
  cnt = last_cell - first_cell + 1 ! cells in this sections

  if ( ElementTypeName(cell_type) .eq. 'QUAD_4' .or. &
       ElementTypeName(cell_type) .eq. 'TRI_3' ) then

    if(parent_flag .eq. 1) then ! If ParentData node in DB is present
      !-------------------------------------------!
      !   Consider boundary conditions            !
      !   with ParentData defined in this block   !
      !-------------------------------------------!
      call Cgns_Mod_Read_2d_Bnd_Section_Connections_With_Parent_Data(  &
        base_id, block_id, sect_id, grid)

      !-----------------------------------------------!
      !   Consider interfaces defined in this block   !
      !-----------------------------------------------!
      call Cgns_Mod_Read_2d_Interface_Section_Connections(  &
        base_id, block_id, sect_id)

    end if ! parent_flag .eq. 1

   end if ! QUAD_4, TRI_3

  !-------------------------------------------------!
  !   Consider three-dimensional cells / sections   !
  !    defined in blocks with unified cell types    !
  !-------------------------------------------------!
  if ( ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PENTA_6') .or.  &
       ( ElementTypeName(cell_type) .eq. 'TETRA_4') ) then

    ! Allocate memory
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_nodes = 8
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_nodes = 5
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_nodes = 6
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_nodes = 4

    allocate(cell_n(n_nodes, cnt))

    call Cg_Elements_Read_F(file_id,       & !(in )
                            base_id,       & !(in )
                            block_id,      & !(in )
                            sect_id,       & !(in )
                            cell_n(1,1),   & !(out)
                            NULL,          & !(out) ! NULL for 3d
                            error)           !(out)

    if (error.ne.0) then
      print "(a)", " # Failed to read section ", sect, " info"
      call Cg_Error_Exit_F()
    endif

    !----------------------------------------------!
    !   Consider boundary conditions               !
    !   with no ParentData defined in this block   !
    !----------------------------------------------!
    ! Search in cell_n(:,:) all QUAD_4/TRI_3 b.c. element connections
    call Cgns_Mod_Read_2d_Bnd_Section_Connections_With_No_Parent_Data( &
      base_id, block_id, sect_id, grid, cell_n)

    ! Fetch received parameters
    do loc = 1, cnt

      ! Calculate absolute cell number
      c = first_cell + loc + cnt_cells - 1

      ! Set the number of nodes for this cell
      grid % cells_n_nodes(c) = n_nodes

      ! Copy individual nodes beloning to this cell
      do n = 1, n_nodes
        grid % cells_n(n, c) = cell_n(n, loc) + cnt_nodes
      end do

      ! Convert from CGNS to Gambit's Neutral file format
      if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) then
        call Swap_Int(grid % cells_n(3, c),  &
                      grid % cells_n(4, c))
        call Swap_Int(grid % cells_n(7, c),  &
                      grid % cells_n(8, c))
      end if
      if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) then
        call Swap_Int(grid % cells_n(3, c),  &
                      grid % cells_n(4, c))
      end if
    end do

    ! Free memory
    deallocate(cell_n)

    if(verbose) then
      print "(a)",   " #================================================="
      print "(a,a26)", " # 3d cell section name: ", trim(sect_name)
      print "(a)",   " #-------------------------------------------------"
      print "(a,i30)", " # Cell section idx: ", sect
      print "(a,a29)", " # Cell section type: ", &
        trim(ElementTypeName(cell_type))
      print "(a,i31)", " # Number of cells: ", cnt
    end if

    if(verbose) then
      print "(a,i26)", " # Corrected first cell: ",  &
        cgns_base(base) % block(block) % section(sect) % first_cell
      print "(a,i27)", " # Corrected last cell: ",  &
        cgns_base(base) % block(block) % section(sect) % last_cell
      print "(a)",     " #-------------------------------------------------"
    end if ! verbose

    ! Count cells in sect
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) cnt_hex = cnt_hex + cnt
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) cnt_pyr = cnt_pyr + cnt
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') cnt_wed = cnt_wed + cnt
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') cnt_tet = cnt_tet + cnt

  end if ! HEXA_8, PYRA_5, PENTA_6, TETRA_4

  !-------------------------------------------------!
  !   Consider three-dimensional cells / sections   !
  !     defined in blocks with mixed cell type      !
  !-------------------------------------------------!
  if(ElementTypeName(cell_type) .eq. 'MIXED') then

    ! Globalized first and last cell of the section
    cgns_base(base) % block(block) % section(sect) % first_cell =  &
    cgns_base(base) % block(block) % section(sect) % first_cell + cnt_cells
    cgns_base(base) % block(block) % section(sect) % last_cell  =  &
    cgns_base(base) % block(block) % section(sect) % last_cell  + cnt_cells

    ! Allocate memory
    allocate(mixed_cell_n(cnt*9))  ! if all are hexa

    call Cg_Elements_Read_F(file_id,           & !(in )
                            base_id,           & !(in )
                            block_id,          & !(in )
                            sect_id,           & !(in )
                            mixed_cell_n(1),   & !(out)
                            NULL,              & !(out) ! NULL for 3d
                            error)               !(out)
    if (error.ne.0) then
      print "(a)", " # Failed to read section ", sect, " info"
      call Cg_Error_Exit_F()
    endif

    ! Fetch received parameters
    loc_hex = 0 ! amount of hex cells in this section
    loc_pyr = 0 ! amount of pyr cells in this section
    loc_wed = 0 ! amount of wed cells in this section
    loc_tet = 0 ! amount of tet cells in this section

    loc = 0
    pos = 0
    do

      ! Increase cell and position counters
      loc = loc + 1
      pos = pos + 1

      ! Exit if end of array has been reached
      if(loc > cnt) exit

      ! Calculate absolute cell number
      c = first_cell + loc + cnt_cells - 1

      ! Increase counters for different types of elements
      cell_type = mixed_cell_n(pos)

      ! Estimate number of nodes for this cell type
      if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) then
        loc_hex = loc_hex + 1
        n_nodes = 8
      end if ! 'HEXA_8'
      if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) then
        loc_pyr = loc_pyr + 1
        n_nodes = 5
      end if ! 'PYRA_5'
      if ( ElementTypeName(cell_type) .eq. 'PENTA_6') then
        loc_wed = loc_wed + 1
        n_nodes = 6
      end if ! 'PENTA_6'
      if ( ElementTypeName(cell_type) .eq. 'TETRA_4') then
        loc_tet = loc_tet + 1
        n_nodes = 4
      end if ! 'TETRA_4'

      ! Store the number of nodes for this cell
      grid % cells_n_nodes(c) = n_nodes

      ! Copy individual nodes beloning to this cell
      do n = 1, n_nodes
        pos = pos + 1
        grid % cells_n(n, c) = mixed_cell_n(pos) + cnt_nodes
      end do ! n = 1, n_nodes

      ! Convert from CGNS to Gambit's Neutral file format
      if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) then
        call Swap_Int(grid % cells_n(3, c),  &
                      grid % cells_n(4, c))
        call Swap_Int(grid % cells_n(7, c),  &
                      grid % cells_n(8, c))
      end if
      if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) then
        call Swap_Int(grid % cells_n(3, c),  &
                      grid % cells_n(4, c))
      end if

    end do ! Fetch received parameters

    !--------------------------------------------------!
    !   Restore cell_n structure                       !
    !     to deal with b.c. in the same way as above   !
    !--------------------------------------------------!
    allocate(cell_n_hex(8, loc_hex))
    allocate(cell_n_pyr(5, loc_pyr))
    allocate(cell_n_wed(6, loc_wed))
    allocate(cell_n_tet(4, loc_tet))

    loc_hex = 0
    loc_pyr = 0
    loc_wed = 0
    loc_tet = 0

    loc = 0
    pos = 0
    do

      ! Increase cell and position counters
      loc = loc + 1
      pos = pos + 1

      ! Exit if end of array has been reached
      if(loc > cnt) exit

      ! Estimate number of nodes for this cell type
      cell_type = mixed_cell_n(pos)
      if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) then
        loc_hex = loc_hex + 1
        cell_n_hex(1:8, loc_hex) = mixed_cell_n(pos+1:pos+8)
        pos = pos + 8
      end if
      if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) then
        loc_pyr = loc_pyr + 1
        cell_n_pyr(1:5, loc_pyr) = mixed_cell_n(pos+1:pos+5)
        pos = pos + 5
      end if
      if ( ElementTypeName(cell_type) .eq. 'PENTA_6') then
        loc_wed = loc_wed + 1
        cell_n_wed(1:6, loc_wed) = mixed_cell_n(pos+1:pos+6)
        pos = pos + 6
      end if
      if ( ElementTypeName(cell_type) .eq. 'TETRA_4') then
        loc_tet = loc_tet + 1
        cell_n_tet(1:4, loc_tet) = mixed_cell_n(pos+1:pos+4)
        pos = pos + 4
      end if
    end do ! mixed_cell_n - > cell_n
    deallocate(mixed_cell_n)

    !--------------------------------------------------------------!
    !   Consider boundary conditions with no ParentData defined    !
    !--------------------------------------------------------------!
    ! Search in cell_n_hex(:,:) all QUAD_4/TRI_3 b.c. element connections
    if (loc_hex > 0) then
      call Cgns_Mod_Read_2d_Bnd_Section_Connections_With_No_Parent_Data( &
        base_id, block_id, sect_id, grid, cell_n_hex)
      deallocate(cell_n_hex)
    end if ! loc_hex > 0

    ! Search in cell_n_pyr(:,:) all QUAD_4/TRI_3 b.c. element connections
    if (loc_pyr > 0) then
      call Cgns_Mod_Read_2d_Bnd_Section_Connections_With_No_Parent_Data( &
        base_id, block_id, sect_id, grid, cell_n_pyr)
      deallocate(cell_n_pyr)
    end if ! loc_pyr > 0

    ! Search in cell_n_wed(:,:) all QUAD_4/TRI_3 b.c. element connections
    if (loc_wed > 0) then
      call Cgns_Mod_Read_2d_Bnd_Section_Connections_With_No_Parent_Data( &
        base_id, block_id, sect_id, grid, cell_n_wed)
      deallocate(cell_n_wed)
    end if ! loc_wed > 0

  ! Search in cell_n_tet(:,:) all QUAD_4/TRI_3 b.c. element connections
    if (loc_tet > 0) then
      call Cgns_Mod_Read_2d_Bnd_Section_Connections_With_No_Parent_Data( &
        base_id, block_id, sect_id, grid, cell_n_tet)
      deallocate(cell_n_tet)
    end if ! loc_tet > 0

    if(verbose) then
      print "(a)",     " #================================================="
      print "(a,a26)", " # 3d cell section name: ", trim(sect_name)
      print "(a)",     " #-------------------------------------------------"
      print "(a,i30)", " # Cell section idx: ", sect
      print "(a,a29)", " # Cell section type: ", &
        trim(ElementTypeName(cell_type))
      print "(a,i31)", " # Number of cells: ", cnt
      print "(a,i26)", " # Corrected first cell: ",  &
        cgns_base(base) % block(block) % section(sect) % first_cell
      print "(a,i27)", " # Corrected last cell: ",  &
        cgns_base(base) % block(block) % section(sect) % last_cell
      print "(a)",     " #-------------------------------------------------"
    end if

    ! Count cells in sect
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) cnt_hex = cnt_hex + loc_hex
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) cnt_pyr = cnt_pyr + loc_pyr
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') cnt_wed = cnt_wed + loc_wed
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') cnt_tet = cnt_tet + loc_tet

  end if ! ElementTypeName(cell_type) .eq. 'MIXED'

  end subroutine
