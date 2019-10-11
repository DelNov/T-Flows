!==============================================================================!
  subroutine Cgns_Mod_Write_Section_Connections(base, block, sect, grid)
!------------------------------------------------------------------------------!
!   Writes elements connection for sect_id [sequential version]                !
!------------------------------------------------------------------------------!
!   Each node in zone/block must have unique id                                !
!   https://cgns.github.io/CGNS_docs_current/midlevel/grid.html                !
!                                                                              !
!   Connections:  | HEXA_8 | PENTA_6 | PYRA_5 | TETRA_4 |                      !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer         :: base, block, sect
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id       ! base index number
  integer           :: block_id      ! block index number
  integer           :: sect_id       ! element section index
  character(len=80) :: sect_name     ! name of the Elements_t node
  integer           :: cell_type     ! types of elements in the section
  integer           :: first_cell    ! index of first element
  integer           :: last_cell     ! index of last element
  integer           :: n_bnd         ! index of last boundary element
  integer           :: error
  integer           :: n_nodes, c, cnt, i, j
  integer           :: cell_n(8, grid % n_nodes)
!==============================================================================!

  ! Set input parameters
  base_id    = base
  block_id   = block
  sect_id    = sect

  sect_name  = trim(cgns_base(base_id)%block(block_id)%section(sect_id)%name)
  cell_type  = cgns_base(base_id)%block(block_id)%section(sect_id)%cell_type
  i          = cgns_base(base_id)%block(block_id)%section(sect_id)%first_cell
  j          = cgns_base(base_id)%block(block_id)%section(sect_id)%last_cell

  ! cells of cell_type
  cnt = j - i + 1
  if ( cnt .ne. 0 ) then

    ! first and last cells have to be shifted according to previous sections
    first_cell = 1 + cnt_cells
    last_cell  = first_cell - 1 + cnt

    n_bnd = 0 ! unsorted boundary elements

    ! Allocate memory
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

    ! Write element data
    call Cg_Section_Write_F(file_id,                   & !(in )
                            base_id,                   & !(in )
                            block_id,                  & !(in )
                            sect_name,                 & !(in )
                            cell_type,                 & !(in )
                            first_cell,                & !(in )
                            last_cell,                 & !(in )
                            n_bnd,                     & !(in )
                            cell_n(1:n_nodes, 1:cnt),  & !(in )
                            sect_id,                   & !(out)
                            error)                       !(out)

    if (error .ne. 0) then
       print*, ' #         Failed to write ', trim(sect_name), ' connections'
       call Cg_Error_Exit_F()
    endif

    ! Print some info
    if(verbose ) then
      print *, '#         ---------------------------------'
      print *, '#         Cell section name: ', sect_name
      print *, '#         ---------------------------------'
      print *, '#         Cell section idx:    ', sect_id
      print *, '#         Cell section type:   ', ElementTypeName(cell_type)
      print *, '#         Number of cells: ', cnt
      print *, '#         First cell:', first_cell
      print *, '#         Last cell: ', last_cell
    end if

    cnt_cells = cnt_cells + cnt
  end if

  end subroutine
