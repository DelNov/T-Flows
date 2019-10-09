!==============================================================================!
  subroutine Cgns_Mod_Read_2d_Interface_Section_Connections(  &
    base_id, block_id, sect_id)
!------------------------------------------------------------------------------!
!   Read 2d section connections and fill interface data
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer              :: base_id, block_id, sect_id
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base       ! base index number
  integer              :: block      ! block index number
  integer              :: sect       ! element section index
  integer              :: int        ! interface index
  character(len=80)    :: sect_name     ! name of the Elements_t node
  character(len=80)    :: int_name      ! name of the interface
  integer              :: int_type      ! type of interface 1-quad, 2-tri, 3-mix
  integer              :: cell_type     ! types of elements in the section
  integer              :: first_cell    ! index of first element
  integer              :: last_cell     ! index of last element
  integer              :: parent_flag
  integer              :: error, int_found 
  integer              :: n_nodes, n, cnt, int_ind, c, j
  integer, allocatable :: face_n(:,:) ! 2d inferface connections
  integer, allocatable :: parent_data(:,:)
!==============================================================================!

  ! Set input parameters
  base  = base_id
  block = block_id
  sect  = sect_id

  ! Introduce some abbreviations
  sect_name   = cgns_base(base) % block(block) % section(sect) % name
  cell_type   = cgns_base(base) % block(block) % section(sect) % cell_type
  first_cell  = cgns_base(base) % block(block) % section(sect) % first_cell
  last_cell   = cgns_base(base) % block(block) % section(sect) % last_cell
  parent_flag = cgns_base(base) % block(block) % section(sect) % parent_flag

  ! Number of cells in this section
  cnt = last_cell - first_cell + 1 ! cells in this sections

  ! Allocate memory
  if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4
  if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3

  allocate(face_n(n_nodes, cnt)); face_n(:,:) = 0
  ! For faces on the boundary of the domain,
  ! the second parent is set to zero
  allocate(parent_data(2*cnt, 2)); parent_data(:,:) = 0

  ! Read element data
  call Cg_Elements_Read_F(file_id,           & !(in )
                          base,              & !(in )
                          block,             & !(in )
                          sect,              & !(in )
                          face_n(:,:),       & !(out)
                          parent_data(:,:),  & !(out)
                          error)               !(out)

  if (error.ne.0) then
    print *, '# Failed to read 2d section connections from: ', sect
    call Cg_Error_Exit_F()
  endif

  do int = 1, cgns_base(base) % block(block) % n_interfaces

    int_found = 0
    do j = 1, cgns_base(base) % block(block) % interface(int) % n_nodes
      if (cgns_base(base) % block(block) % interface(int) % &
          belongs_to_sect(j) .eq. sect) then
        int_found = int_found + 1
      end if ! if face is in sect
    end do ! j = 1, cgns_base(base) % block(block) % interface(int) % n_nodes

    if (int_found .ne. 0) then

      int_name = trim(cgns_base(base) % block(block) % interface(int) % name)
      int_ind = cgns_base(base) % block(block) % interface(int) % id
      int_type = cgns_base(base) % block(block) % interface(int) % int_type

      if(verbose) then
        print '(a)',     ' #-----------------------------------------------'
        print '(a,a24)', ' # 2d cell section name: ', trim(sect_name)
        print '(a)',     ' #-----------------------------------------------'
        print '(a,a30)', ' # Interface name: ', trim(int_name)
        print '(a,i29)', ' # Interface index: ', &
          cgns_base(base) % block(block) % interface(int) % id
        print '(a,i31)', ' # Section index: ', sect
        print '(a,a30)', ' # Interface type: ', &
          trim(ElementTypeName(cell_type))
        print '(a,l25)', ' # Marked for deletion: ', cgns_base(base) % &
          block(block) % interface(int) % marked_for_deletion
        !print '(a)',     ' # Interface cells connection table (sample): '
        !do loc = 1, min(6, cnt)
        !  print '(a,a16,4i8)', ' # ', ' ', &
        !    (face_n(n,loc), n = 1, n_nodes)
        !end do
        !print *, '#       Interface parent data (sample): '
        !do loc = 1, min(6, cnt)
        !    print '(a9,8i7)', ' ',parent_data(loc, 1)
        !end do
      end if ! verbose

      ! If interface is not marked for deletion -> unique
      if ( .not. cgns_base(base) % block(block) % interface(int) % &
        marked_for_deletion) then

        ! Add unique interface (considering mixed)
        if (int_type <    3) cnt_int = cnt_int + 2
        if (int_type .eq. 3) cnt_int = cnt_int + 1

        do c = 1, cnt
          do n = 1, n_nodes
            interface_cells(1, c + cnt_int_cells, n, int_ind) = &
              face_n(n, c) + cnt_nodes
          end do ! n = 1, n_nodes
        end do ! c = 1, cnt

      else ! marked for deletion

        do c = 1, cnt
          do n = 1, n_nodes
            interface_cells(2, c + cnt_int_cells, n, int_ind) = &
              face_n(n, c) + cnt_nodes
          end do ! n = 1, n_nodes
        end do ! c = 1, cnt

      end if ! .not. marked for deletion

      ! Move starting point for next interface
      cnt_int_cells = cnt_int_cells + cnt
    end if ! int_found
  end do ! int = 1, cgns_base(base) % block(block) % n_interfaces

  deallocate(face_n)
  deallocate(parent_data)

  end subroutine
