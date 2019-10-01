!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Connections(base, block, sect, grid,  &
                                               parent_flag)
!------------------------------------------------------------------------------!
!   Read elements connection info for current sect
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer         :: base, block, sect
  integer         :: parent_flag
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base_id       ! base index number
  integer              :: block_id      ! block index number
  integer              :: sect_id       ! element section index
  character(len=80)    :: sect_name     ! name of the Elements_t node
  character(len=80)    :: int_name      ! name of the interface
  character(len=80)    :: bnd_name      ! name of the interface
  integer              :: min_name_l
  integer              :: int_type      ! type of interface 1-quad, 2-tri, 3-mix
  integer              :: cell_type     ! types of elements in the section
  integer              :: first_cell    ! index of first element
  integer              :: last_cell     ! index of last element
  integer              :: error
  integer              :: n_nodes, loc, c, n, cell, dir, cnt, bc, int,  &
                          int_id, pos
  integer, allocatable :: cell_n(:,:), mixed_cell_n(:)
  integer, allocatable :: face_n(:,:)
  integer, allocatable :: interface_n(:,:)
  integer, allocatable :: parent_data(:,:)
  integer              :: parent_datum = 0  ! for cells there are no parents
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

  ! For faces on the boundary of the domain, the second parent is set to zero"
  allocate(parent_data(2*cnt, 2))

  !--------------------------------------------------------!
  !   Consider boundary conditions defined in this block   !
  !--------------------------------------------------------!
  do bc = 1, cgns_base(base) % block(block) % n_bnd_conds

    bnd_name    = trim(cgns_base(base) % block(block) % bnd_cond(bc) % name)
    min_name_l  = min(len(trim(bnd_name)), len(trim(sect_name)))

    if(bnd_name(1:min_name_l) .eq. sect_name(1:min_name_l)) then

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
      !if ( ElementTypeName(cell_type) .eq. 'QUAD_4') cnt_qua = cnt_qua + cnt
      !if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) cnt_tri = cnt_tri + cnt

      ! Allocate memory
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3
      allocate(face_n(n_nodes, cnt))

      call Cg_Elements_Read_F(file_id,      & !(in )
                              base_id,      & !(in )
                              block_id,     & !(in )
                              sect_id,      & !(in )
                              face_n(1,1),  & !(out)
                              parent_data,  & !(out)
                              error)          !(out)

      ! Fetch the data if parent is provided
      if(parent_flag .eq. 1) then
        do loc = 1, cnt
          cell = parent_data(loc, 1) + cnt_cells
          dir  = parent_data(loc, 2)
          grid % cells_bnd_color(dir,cell) =  &
               cgns_base(base) % block(block) % bnd_cond(bc) % color
        end do

      ! Parent data not provided
      else
        do loc = 1, cnt
          grid % cells_n_nodes(-cnt_bnd_cells-loc) = n_nodes

          ! Copy individual nodes beloning to this cell
          do n = 1, n_nodes
            grid % cells_n(n, -cnt_bnd_cells-loc) =  &
                    face_n(n, loc) + cnt_nodes
          end do

          grid % bnd_cond % color(-cnt_bnd_cells-loc) =  &
                    cgns_base(base) % block(block) % bnd_cond(bc) % color
        end do
      end if

      ! Update number of boundary cells in the block
      cnt_bnd_cells = cnt_bnd_cells + cnt

      if(verbose) then
        print *, "#         Connection table (sample): "
        do loc = 1, min(6, cnt)
          print '(a8,4i7)', " ", (face_n(n,loc), n = 1, n_nodes)
        end do
        if(parent_flag .eq. 1) then
          print *, "#         Parent data (sample): "
          do loc = 1, min(6, cnt)
            print '(a10,2i7)', " ", parent_data(loc, 1), parent_data(loc, 2)
          end do
        end if
      end if

      deallocate(face_n)

    end if
  end do

  !-----------------------------------------------!
  !   Consider interfaces defined in this block   !
  !-----------------------------------------------!
  do int = 1, cgns_base(base) % block(block) % n_interfaces

    int_name = trim(cgns_base(base) % block(block) % interface(int) % name)
    int_id = cgns_base(base) % block(block) % interface(int) % id
    int_type = cgns_base(base) % block(block) % interface(int) % int_type

    if(index(trim(sect_name), trim(int_name), back = .true.) .ne. 0) then

      ! Allocate memory
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3
      allocate(interface_n(n_nodes, cnt))  

      call Cg_Elements_Read_F(file_id,           & !(in )
                              base_id,           & !(in )
                              block_id,          & !(in )
                              sect_id,           & !(in )
                              interface_n(1,1),  & !(out)
                              parent_data,       & !(out)
                              error)               !(out)

      ! If interface is not marked for deletion
      if ( .not. cgns_base(base) % &
        block(block) % interface(int) % marked_for_deletion) then

        ! Add unique interface (considering mixed)
        if (int_type <    3) cnt_int = cnt_int + 2
        if (int_type .eq. 3) cnt_int = cnt_int + 1

        ! Fetch first interface
        do loc = 1, cnt
          do n = 1, n_nodes
          interface_cells(1, loc + cnt_int_cells, n, int_id) = &
            interface_n(n, loc) + cnt_nodes
          end do ! n
        end do ! loc

      else

        ! Fetch second interface
        do loc = 1, cnt
          do n = 1, n_nodes
          interface_cells(2, loc + cnt_int_cells, n, int_id) = &
            interface_n(n, loc) + cnt_nodes
          end do ! n
        end do ! loc

      end if

      ! Fetch received parameters
      cnt_int_cells = cnt_int_cells + cnt

      if(verbose) then
        print *, '#         ---------------------------------'
        print *, '#         Interface name:  ', trim(sect_name)
        print *, '#         ---------------------------------'
        print *, '#         Interface index: ', cgns_base(base) % &
        block(block) % interface(int) % id
        print *, '#         Section index: ', sect
        print *, '#         Interface type:  ', ElementTypeName(cell_type)
        print *, '#         Marked for deletion:  ', cgns_base(base) % &
        block(block) % interface(int) % marked_for_deletion
        print *, "#         Interface cells connection table (sample): "
        do loc = 1, min(6, cnt)
          print '(a9,8i7)', " ", (interface_n(n,loc), n = 1, n_nodes)
        end do
        print *, "#         Interface parent data (sample): "
        do loc = 1, min(6, cnt)
            print '(a9,8i7)', " ",parent_data(loc, 1)
        end do
      end if

        deallocate(interface_n)
    end if
  end do

  !-------------------------------------------------!
  !   Consider three-dimensional cells / sections   !
  !    defined in blocks with unified cell types    !
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

    call Cg_Elements_Read_F(file_id,       & !(in )
                            base_id,       & !(in )
                            block_id,      & !(in )
                            sect_id,       & !(in )
                            cell_n(1,1),   & !(out)
                            parent_datum,  & !(out)
                            error)           !(out)

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
        print *, "#         Connection table (sample): "
      do loc = 1, min(6, cnt)
        print '(a9,8i7)', " ", (cell_n(n,loc), n = 1, n_nodes)
      end do
    end if

    deallocate(cell_n)
  end if

  !-------------------------------------------------!
  !   Consider three-dimensional cells / sections   !
  !     defined in blocks with mixed cell types     !
  !-------------------------------------------------!
  if(ElementTypeName(cell_type) .eq. 'MIXED') then

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

    ! Allocate memory
    allocate(mixed_cell_n(cnt*9))  ! if all are hexa

    call Cg_Elements_Read_F(file_id,           & !(in )
                            base_id,           & !(in )
                            block_id,          & !(in )
                            sect_id,           & !(in )
                            mixed_cell_n(1),   & !(out)
                            parent_datum,      & !(out)
                            error)               !(out)

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
      if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) cnt_hex = cnt_hex + cnt
      if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) cnt_pyr = cnt_pyr + cnt
      if ( ElementTypeName(cell_type) .eq. 'PENTA_6') cnt_wed = cnt_wed + cnt
      if ( ElementTypeName(cell_type) .eq. 'TETRA_4') cnt_tet = cnt_tet + cnt

      ! Estimate number of nodes for this cell type
      cell_type = mixed_cell_n(pos)
      if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_nodes = 8
      if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_nodes = 5
      if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_nodes = 6
      if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_nodes = 4

      ! Store the number of nodes for this cell
      grid % cells_n_nodes(c) = n_nodes

      ! Copy individual nodes beloning to this cell
      do n = 1, n_nodes
        pos = pos + 1
        grid % cells_n(n, c) = mixed_cell_n(pos) + cnt_nodes
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

    deallocate(mixed_cell_n)
  end if

  end subroutine
