!==============================================================================!
  subroutine Cgns_Mod_Write_Section_Connections(base, block, sect, grid)
!------------------------------------------------------------------------------!
!   Writes elements connection for sect_id [parallel version]                  !
!------------------------------------------------------------------------------!
!   Each node in zone/block must have unique id
!   https://cgns.github.io/CGNS_docs_current/midlevel/grid.html                !
!------------------------------------------------------------------------------!
!   Array structures in current function are strictly followings:              !
!                                                                              !
!   Cell type:    |      HEXA_8      |     PENTA_6      |       PYRA_5     |...!
!   Connections:  |-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|...!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer         :: base, block, sect
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: base_id       ! base index number
  integer       :: block_id      ! block index number
  integer       :: sect_id       ! element section index
  character(SL) :: sect_name     ! name of the Elements_t node
  integer       :: cell_type     ! types of elements in the section
  integer       :: n_bnd         ! index of last boundary element
  integer       :: error
  integer       :: n_nodes, c, cnt, i, j
  integer       :: first_cell    ! look at array structure at the header
  integer       :: last_cell     ! look at array structure at the header
  integer       :: cell_n(8, grid % n_nodes)
!==============================================================================!

  ! Set input parameters
  base_id    = base
  block_id   = block
  sect_id    = sect

  sect_name  = trim(cgns_base(base_id)%block(block_id)%section(sect_id)%name)
  cell_type  = cgns_base(base_id)%block(block_id)%section(sect_id)%cell_type
  i          = cgns_base(base_id)%block(block_id)%section(sect_id)%first_cell
  j          = cgns_base(base_id)%block(block_id)%section(sect_id)%last_cell

  ! cells of cell_type on this_proc
  cnt = j - i + 1
  if ( cnt .ne. 0 ) then
  !-------------------------------------------!
  !   Create empty sect_name node in DB block !
  !-------------------------------------------!

    ! total cells of cell_type
    c = cnt
    call Comm_Mod_Global_Sum_Int(c)
    n_bnd = 0 ! unsorted boundary elements

    ! Create empty elem node in DB
    call Cgp_Section_Write_F(file_id,        & !(in )
                             base_id,        & !(in )
                             block_id,       & !(in )
                             sect_name,      & !(in )
                             cell_type,      & !(in )
                             1 + cnt_cells,  & !(in )
                             c + cnt_cells,  & !(in )
                             n_bnd,          & !(in )
                             sect_id,        & !(out)
                             error)            !(out)
    if (error .ne. 0) then
      print*, '*FAILED* to create empty ', trim(sect_name)
      call Cgp_Error_Exit_F()
    endif

    ! Print some info
    if(verbose .and. this_proc.lt.2) then
      print *, '#         ---------------------------------'
      print *, '#         Cell section name: ', sect_name
      print *, '#         ---------------------------------'
      print *, '#         Cell section idx:    ', sect_id
      print *, '#         Cell section type:   ', ElementTypeName(cell_type)
      print *, '#         Cell section first cell:   ', 1 + cnt_cells
      print *, '#         Cell section last cell:   ',  c + cnt_cells
    end if


  !------------------------------------------------!
  !   Mapping 1:nc -> Connection structure above   !
  !------------------------------------------------!

    i = cnt ! cnt_hex/pyr/wed/tet on this_proc
    call Cgns_Mod_Get_Arrays_Dimensions(j, i)

    first_cell = j + cnt_cells
    last_cell  = first_cell - 1 + cnt

    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_nodes = 8
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_nodes = 6
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_nodes = 5
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_nodes = 4

    ! Convert T-FlowS -> CGNS [same as VTK]
    i = 1
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      if (grid % cells_n_nodes(c).eq.8 .and. n_nodes.eq.8 ) then ! hex
        cell_n (1, i) = grid % cells_n(1, c)
        cell_n (2, i) = grid % cells_n(2, c)
        cell_n (3, i) = grid % cells_n(4, c)
        cell_n (4, i) = grid % cells_n(3, c)
        cell_n (5, i) = grid % cells_n(5, c)
        cell_n (6, i) = grid % cells_n(6, c)
        cell_n (7, i) = grid % cells_n(8, c)
        cell_n (8, i) = grid % cells_n(7, c)
        i = i + 1
      elseif (grid % cells_n_nodes(c).eq.6 .and. n_nodes.eq.6) then ! wedge
        cell_n (1:6, i) = grid % cells_n(1:6, c)
        i = i + 1
      elseif (grid % cells_n_nodes(c).eq.5 .and. n_nodes.eq.5) then ! pyramid
        cell_n (1, i) = grid % cells_n(5, c)
        cell_n (2, i) = grid % cells_n(1, c)
        cell_n (3, i) = grid % cells_n(2, c)
        cell_n (4, i) = grid % cells_n(4, c)
        cell_n (5, i) = grid % cells_n(3, c)
        i = i + 1
      elseif (grid % cells_n_nodes(c).eq.4 .and. n_nodes.eq.4) then ! tetra
        cell_n (1:4, i) = grid % cells_n(1:4, c)
        i = i + 1
      end if
    end do

    ! Shift cell_n value according to a number of nodes in coord array
    i = grid % n_nodes
    call Cgns_Mod_Get_Arrays_Dimensions(j, i)
    cell_n = cell_n + j - 1

    !--------------------------------------!
    !   Fill empty sect_name in DB block   !
    !--------------------------------------!

    ! Fill that node with grid coordinates
    call Cgp_Elements_Write_Data_F(file_id,                   & !(in )
                                   base_id,                   & !(in )
                                   block_id,                  & !(in )
                                   sect_id,                   & !(in )
                                   first_cell,                & !(in )
                                   last_cell,                 & !(in )
                                   cell_n(1:n_nodes, 1:cnt),  & !(in )
                                   error)                       !(out)

    if (error .ne. 0) then
      print *, '*FAILED* to fill ', sect_id
      call Cgp_Error_Exit_F()
    endif

    ! Print some info
    if(verbose) then
      print *, '#         First cell:', first_cell, ' (P:',this_proc,')'
      print *, '#         Last cell: ', last_cell,  ' (P:',this_proc,')'
    end if

    c = cnt
    call Comm_Mod_Global_Sum_Int(c)
    cnt_cells = cnt_cells + c

  end if

  end subroutine
