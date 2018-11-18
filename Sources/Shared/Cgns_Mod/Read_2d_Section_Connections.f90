!==============================================================================!
  subroutine Cgns_Mod_Read_2d_Section_Connections(base, block, sect, grid)
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
  character(len=80)    :: int_name      ! name of the interface
  integer              :: int_type      ! type of interface 1-quad, 2-tri, 3-mix
  integer              :: cell_type     ! types of elements in the section
  integer              :: first_cell    ! index of first element
  integer              :: last_cell     ! index of last element
  integer              :: parent_flag
  integer              :: error
  integer              :: n_nodes, loc, c, n, cell, dir, cnt, bc, int, int_id
  integer              :: i, j, k
  integer, allocatable :: cell_n(:,:)
  !integer(kind=4), allocatable :: face_n(:,:)
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

  if ( ElementTypeName(cell_type) .eq. 'QUAD_4' .or. &
       ElementTypeName(cell_type) .eq. 'TRI_3' ) then
    !--------------------------------------------------------!
    !   Consider boundary conditions defined in this block   !
    !--------------------------------------------------------!
    do bc = 1, cgns_base(base) % block(block) % n_bnd_conds

      ! Allocate memory
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3

      ! If ParentElements node is present
      if(parent_flag .eq. 1) then

      !"For faces on the boundary of the domain, the second parent is set to zero"
        allocate(parent_data(2*cnt, 2)) ! provided in Read_Section_Info
        allocate(face_n(n_nodes, cnt))

        ! Read element data
        call Cg_Elements_Read_F(file_id,      & !(in )
                                base_id,      & !(in )
                                block_id,     & !(in )
                                sect_id,      & !(in )
                                face_n,  & !(out)
                                parent_data,  & !parent data
                                error)          !(out)

        ! Fetch the data
        do loc = 1, cnt
          cell = parent_data(loc, 1) + cnt_cells
          dir  = parent_data(loc, 2)
          grid % cells_bnd_color(dir, cell) =  &
               cgns_base(base) % block(block) % bnd_cond(bc) % color
        end do

        if(verbose) then
          print *, "#         Connection table (sample): "
          do loc = 1, min(6, cnt)
            print '(a8,4i7)', " ", (face_n(n,loc), n = 1, n_nodes)
          end do
          print *, "#         Parent data (sample): "
          do loc = 1, min(6, cnt)
            print '(a10,2i7)', " ", parent_data(loc, 1), parent_data(loc, 2)
          end do
        end if

        deallocate(face_n)
      else

        allocate(face_n(n_nodes, cnt))

        ! Read element data
        call Cg_Elements_Read_F(file_id,  & !(in )
                                base_id,  & !(in )
                                block_id, & !(in )
                                sect_id,  & !(in )
                                face_n,   & !(out)
                                NULL,     & !parent data
                                error)      !(out)

      ! If point of b.c. is inside this section -> this section is b.c.
      k = 0
      do j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes
        ! Fetch received parameters
        do loc = 1, cnt
          if (loc .eq. &
            cgns_base(base) % block(block) % bnd_cond(bc) % point_list(j)) then

            ! Set the number of nodes for this cell
            cgns_base(base) % block(block) % bnd_cond(bc) % cells_n_nodes(j) = &
              n_nodes

            ! Copy individual nodes beloning to this cell
            do n = 1, n_nodes
              cgns_base(base) % block(block) % bnd_cond(bc) % cells_n_nodes(j) = face_n(n, loc) + cnt_nodes
          end if
        end do
      end do

      if(verbose .and. k .ne. 0) then
        print *, '#         ---------------------------------'
        print *, '#         Bnd section name:  ', trim(sect_name)
        print *, '#         ---------------------------------'
        print *, '#         Bnd section index: ', sect
        print *, '#         Bnd section type: ', ElementTypeName(cell_type)
        print *, '#         Bnd condition color: ',   &
                 cgns_base(base) % block(block) % bnd_cond(bc) % color
        print *, '#         Number of faces: ', k
        print *, '#         They belong to b.c.: ', &
          trim(cgns_base(base) % block(block) % bnd_cond(bc) % name)
      end if

      print *, "n_nodes", n_nodes
      print *, "cnt", cnt
      print *, "face_n:"
      do loc = 1, min(6, cnt)
        print '(a8,4i7)', " ", (face_n(n,loc), n = 1, n_nodes)
      end do

      ! Update number of boundary cells in the block
      cnt_block_bnd_cells = cnt_block_bnd_cells + k

        ! Fetch the data
        do loc = 1, cnt
          cell = parent_data(loc, 1) + cnt_cells
          dir  = parent_data(loc, 2)
          grid % cells_bnd_color(dir, cell) =  &
               cgns_base(base) % block(block) % bnd_cond(bc) % color
        end do

        deallocate(face_n)
        stop
      end if ! parent_flag .eq. 1

    end do ! bc

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

        ! Read element data
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
  end if ! QUAD_4 and TRI_3

  end subroutine