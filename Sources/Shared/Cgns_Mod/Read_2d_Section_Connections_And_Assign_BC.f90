!==============================================================================!
  subroutine Cgns_Mod_Read_2d_Section_Connections_And_Assign_BC(  &
    base, block, sect, grid, cell_n)
!------------------------------------------------------------------------------!
!   Read 2d elements connection info for current sect
!   Apply b.c. and handle interfaces be searching 2d cell inside 3d cell
!   Prerequisites: all 3d cells in this block are already in grid structure
!   https://cgns.github.io/CGNS_docs_current/sids/conv.html
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer             :: base, block, sect
  integer, intent(in) :: cell_n(:,:)
  type(Grid_Type)     :: grid
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
  !!integer(kind=4), allocatable :: face_n(:,:)
  integer, allocatable :: face_n(:,:)
  integer, allocatable :: interface_n(:,:)
  integer, allocatable :: parent_data(:,:)
  integer              :: face_quad(4)
  integer              :: face_tri(3)
  ! cgns HEXA_8 cell faces nodal connections
  integer, parameter, dimension(6, 4) :: face_hexa_8 =                                     reshape( (/ 1, 4, 3, 2,  &
                                                     1, 2, 6, 5,  &
                                                     2, 3, 7, 6,  &
                                                     3, 4, 8, 7,  &
                                                     1, 5, 8, 4,  &
                                                     5, 6, 7, 8  /), (/6, 4/) )
  ! cgns PYRA_5 cell faces nodal connections
  integer, parameter, dimension(6, 4) :: face_pyra_5 =                                     reshape( (/ 1, 4, 3, 2,  &
                                                     1, 2, 5,-1,  &
                                                     2, 3, 5,-1,  &
                                                     3, 4, 5,-1,  &
                                                     4, 1, 5,-1,  &
                                                    -1,-1,-1,-1  /), (/6, 4/) )
  ! cgns PENTA_6 cell faces nodal connections
  integer, parameter, dimension(6, 4) :: face_penta_6 =                                     reshape( (/ 1, 2, 5, 4,  &
                                                     2, 3, 6, 5,  &
                                                     3, 1, 4, 6,  &
                                                     1, 3, 2,-1,  &
                                                     4, 5, 6,-1,  &
                                                    -1,-1,-1,-1  /), (/6, 4/) )
  ! cgns TETRA_4 cell faces nodal connections
  integer, parameter, dimension(6, 4) :: face_tetra_4 =                                     reshape( (/ 1, 3, 2,-1,  &
                                                     1, 2, 4,-1,  &
                                                     2, 3, 4,-1,  &
                                                     3, 1, 4,-1,  &
                                                    -1,-1,-1,-1,  &
                                                    -1,-1,-1,-1  /), (/6, 4/) )
  !integer :: face (6,4)
  integer, allocatable :: bc_face(:)
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

        print *, "trim(sect_name) ",trim(sect_name)
        print *, "is inside"
        print *, "trim(cgns_base(base) % block(block) % bnd_cond(bc) % name ", &
        trim(cgns_base(base) % block(block) % bnd_cond(bc) % name)

      ! Allocate memory
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3

      ! Array to store 2d cells connections
      allocate(face_n(n_nodes, cnt))

      ! If ParentElements node is present (provided in Read_Section_Info)
      if(parent_flag .eq. 1) then

        !"For faces on the boundary of the domain,
        ! the second parent is set to zero"
        allocate(parent_data(2*cnt, 2))

        ! Read element data
        call Cg_Elements_Read_F(file_id,      & !(in )
                                base_id,      & !(in )
                                block_id,     & !(in )
                                sect_id,      & !(in )
                                face_n,       & !(out)
                                parent_data,  & !(out)
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

     ! no Parent Data is present
      else

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
        do i = first_cell, last_cell
          do j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes

            if( i .eq. &
              cgns_base(base)%block(block) % bnd_cond(bc) % point_list(j) ) then
              ! this b.c. cells belongs to section:
            end if

          end do
        end do

        allocate(bc_face(n_nodes))

        do i = 1, n_nodes ! from 1 to n_faces
         bc_face(:) = face_n(:, i)
          do j = 1, size(cell_n, 2) ! from 1 to n_cells
            if     ( size(cell_n, 1) .eq. 8 ) then ! HEXA_8
              print *, '! hex: '
            elseif ( size(cell_n, 1) .eq. 5 ) then ! PYRA_5
              print *, '! pyra: '
            elseif ( size(cell_n, 1) .eq. 6 ) then ! PENTA_6
              print *, '! penta: '
            elseif ( size(cell_n, 1) .eq. 4 ) then ! TETRA_4
              print *, '! tetra: '
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

        print *, "face_n:"
        do loc = 1, min(6, cnt)
          print '(a8,4i7)', " ", (face_n(n,loc), n = 1, n_nodes)

        end do

        print *, "cell_n:"
        do loc = 1, min(6, cnt)
          print '(a9,8i7)', " ", (cell_n(n,loc), n = 1, size(cell_n, 1))

        end do

        stop


        if(verbose .and. k .ne. 0) then
          !do loc = 1, grid % cells_n_nodes(c)
            !grid % cells_n(n, c)
          !end do
          
        end if

        deallocate(face_n)

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