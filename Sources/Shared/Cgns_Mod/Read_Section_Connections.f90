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
  integer              :: n_nodes, loc, c, n, cnt
  integer, allocatable :: cell_n(:,:)
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

    !------------------------------------------------------------------------!
    !   Consider boundary conditions with ParentData defined in this block   !
    !------------------------------------------------------------------------!

    if(parent_flag .eq. 1) then ! If ParentData node in DB is present
      call Cgns_Mod_Read_2d_Bnd_Section_Connections_With_Parent_Data(  &
        base_id, block_id, sect_id, grid)

    !-----------------------------------------------!
    !   Consider interfaces defined in this block   !
    !-----------------------------------------------!
      call Cgns_Mod_Read_2d_Interface_Section_Connections(  &
        base_id, block_id, sect_id)
    end if ! if parent_flag .eq. 1

   end if ! QUAD_4, TRI_3

  !-------------------------------------------------!
  !   Consider three-dimensional cells / sections   !
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
                            cell_n(:,:),   & !(out)
                            NULL,          & !(out) ! NULL for 3d
                            error)           !(out)

    if (error.ne.0) then
      print "(a)", " # Failed to read section ", sect, " info"
      call Cg_Error_Exit_F()
    endif
    !--------------------------------------------------------------!
    !   Consider boundary conditions with no ParentData defined    !
    !--------------------------------------------------------------!
    ! Search in cell_n(:,:) all QUAD_4/TRI_3 b.c. element connections
    call Cgns_Mod_Read_2d_Bnd_Section_Connections_With_No_Parent_Data( &
      base_id, block_id, sect_id, grid, cell_n)


    ! Fetch received parameters
    do loc = 1, cnt

      ! Calculate absolute cell number
      c = loc - 1 + first_cell

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
      print "(a)", " #-------------------------------------------------"
      print "(a,a26)", ' # 3d cell section name: ', trim(sect_name)
      print "(a)", " #-------------------------------------------------"
      print "(a,i30)", ' # Cell section idx: ', sect
      print "(a,a29)", ' # Cell section type: ', &
        trim(ElementTypeName(cell_type))
      print "(a,i31)", ' # Number of cells: ', cnt
    end if

    if(verbose) then
      print "(a,i26)", " # Corrected first cell: ",  &
        cgns_base(base) % block(block) % section(sect) % first_cell
      print "(a,i27)", " # Corrected last cell: ",  &
        cgns_base(base) % block(block) % section(sect) % last_cell
      print "(a)", " #-------------------------------------------------"
    end if ! verbose

    ! Count cells in sect
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) cnt_hex = cnt_hex + cnt
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) cnt_pyr = cnt_pyr + cnt
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') cnt_wed = cnt_wed + cnt
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') cnt_tet = cnt_tet + cnt

  end if ! HEXA_8, PYRA_5, PENTA_6, TETRA_4

  end subroutine