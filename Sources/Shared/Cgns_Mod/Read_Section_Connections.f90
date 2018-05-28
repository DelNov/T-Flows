!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Connections(base, block, sect, grid)
!------------------------------------------------------------------------------!
!   Read elements connection info for current sect
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
  integer              :: n_nodes, loc, c, n, cell, dir, cnt, bc
  integer, allocatable :: cell_n(:,:)
  integer, allocatable :: face_n(:,:)
  integer, allocatable :: parent_data(:,:)
  integer              :: parent_datum = 0  ! for cells there are no parents
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  sect_id  = sect

  ! Introduce some abbraviations 
  sect_name   = cgns_base(base) % block(block) % section(sect) % name
  cell_type   = cgns_base(base) % block(block) % section(sect) % cell_type
  first_cell  = cgns_base(base) % block(block) % section(sect) % first_cell
  last_cell   = cgns_base(base) % block(block) % section(sect) % last_cell
  parent_flag = cgns_base(base) % block(block) % section(sect) % parent_flag

  ! Number of cells in this section
  cnt = last_cell - first_cell + 1 ! cells in this sections

  if(parent_flag .eq. 1) then
    allocate(parent_data(2*cnt,2))
  end if

  !--------------------------------------------------------!
  !   Consider boundary conditions defined in this block   !
  !--------------------------------------------------------!
  do bc = 1, cgns_base(base) % block(block) % n_bnd_conds
    if(sect_name .eq. cgns_base(base) % block(block) % bnd_cond(bc) % name) then

      if(verbose) then
        print *, '#         ---------------------------------'
        print *, '#         Bnd section name:  ', sect_name
        print *, '#         ---------------------------------'
        print *, '#         Bnd section index:  ', sect
        print *, '#         Bnd section type:   ', ElementTypeName(cell_type)
        print *, '#         Bnd condition color:',   &
                 cgns_base(base) % block(block) % bnd_cond(bc) % color
        print *, '#         Number of faces:    ', cnt
      end if

      ! Count boundary cells
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') cnt_qua = cnt_qua + cnt
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) cnt_tri = cnt_tri + cnt

      ! Update numer of boundary cells in the block
      cnt_block_bnd_cells = cnt_block_bnd_cells + cnt

      ! Allocate memory
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3
      allocate(face_n(n_nodes, cnt))

      call Cg_Elements_Read_F(file_id,       & 
                              base_id,       &
                              block_id,      &
                              sect_id,       &
                              face_n(1,1),   &
                              parent_data,   &
                              error)          
      ! Fetch the data
      do loc=1,cnt  ! I have no clue why the size has to be 2*cnt
        cell = parent_data(loc,1) + cnt_cells
        dir  = parent_data(loc,2)
        grid % cells_bnd_color(dir,cell) =  &
             cgns_base(base) % block(block) % bnd_cond(bc) % color
      end do
     
      if(verbose) then
        do loc = 1,min(8,cnt)
          print '(4i7)', (face_n(n,loc), n = 1, n_nodes)
        end do
        do loc = 1,min(8,cnt)
          print '(a2,3i7)', 'loc=', loc, parent_data(loc,1), parent_data(loc,2)
        end do
      end if

      deallocate(face_n)

    end if
  end do

  !-------------------------------------------------!
  !   Consider three-dimensional cells / sections   !
  !-------------------------------------------------!
  if ( ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PENTA_6') .or.  &
       ( ElementTypeName(cell_type) .eq. 'TETRA_4') ) then

    if(verbose) then
      print *, '#         ---------------------------------'
      print *, '#         Cell section name: ', sect_name
      print *, '#         ---------------------------------'
      print *, '#         Cell section idx:    ', sect
      print *, '#         Cell section type:   ', ElementTypeName(cell_type)
      print *, '#         Number of cells:     ', cnt
      print *, '#         Corrected first cell:', first_cell + cnt_cells
      print *, '#         Corrected last cell: ', last_cell  + cnt_cells
    end if

    ! Globalized first and last cell of the section
    cgns_base(base) % block(block) % section(sect) % first_cell =  &
    cgns_base(base) % block(block) % section(sect) % first_cell + cnt_cells
    cgns_base(base) % block(block) % section(sect) % last_cell  =  &
    cgns_base(base) % block(block) % section(sect) % last_cell  + cnt_cells

    ! Count cells in sect
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) cnt_hex = cnt_hex + cnt
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) cnt_pyr = cnt_pyr + cnt
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') cnt_wed = cnt_wed + cnt
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') cnt_tet = cnt_tet + cnt

    ! Allocate memory
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_nodes = 8
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_nodes = 5
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_nodes = 6
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_nodes = 4
    allocate(cell_n(n_nodes, cnt))

    call Cg_Elements_Read_F(file_id,       & 
                            base_id,       &
                            block_id,      &
                            sect_id,       &
                            cell_n(1,1),   &
                            parent_datum,  &
                            error)          
    
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

    if(verbose) then
      do loc = 1, min(8,cnt)
        print '(8i7)', (cell_n(n,loc), n = 1, n_nodes)
      end do
    end if
 
    deallocate(cell_n)
  end if

  end subroutine
