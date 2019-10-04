!==============================================================================!
  subroutine Cgns_Mod_Read_2d_Bnd_Section_Connections_With_Parent_Data(  &
    base_id, block_id, sect_id, grid)
!------------------------------------------------------------------------------!
!   Read 2d elements connections for current 3d sect                           !
!   ParentData node in DB is assumed to be present                             !
!   It is used to assign b.c.                                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer              :: base_id, block_id, sect_id
  type(Grid_Type)      :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base          ! base index number
  integer              :: block         ! block index number
  integer              :: sect          ! element section index
  integer              :: bc            ! boundary color id
  character(len=80)    :: sect_name     ! name of the Elements_t node
  integer              :: cell_type     ! types of elements in the section
  integer              :: first_cell    ! index of first element
  integer              :: last_cell     ! index of last element
  integer              :: error
  integer              :: n_nodes, loc, n, cnt, j
  integer, allocatable :: parent_data(:,:)
  integer, allocatable :: face_n(:,:)
  integer              :: dir, cell, bc_found
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

  ! Number of cells in this section
  cnt = last_cell - first_cell + 1 ! cells in this sections

  ! Allocate memory
  if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4
  if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3

  allocate(face_n(n_nodes, cnt)); face_n(:,:) = 0
  !"For faces on the boundary of the domain,
  ! the second parent is set to zero"
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
    print "(a)", " # Failed to read 2d section connections from: ", sect
    call Cg_Error_Exit_F()
  endif

  !----------------------!
  !   Count b.c. faces   !
  !----------------------!

  do bc = 1, cgns_base(base) % block(block) % n_bnd_conds

    bc_found = 0
    do j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes
      if (cgns_base(base) % block(block) % bnd_cond(bc) % &
        belongs_to_sect(j) .eq. sect_id) then
        bc_found = bc_found + 1
      end if
    end do ! j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes

    if(bc_found .ne. 0) then
      if(verbose) then
        print "(a)",     " #==============================================="
        print "(a,a24)", " # 2d cell section name: ", trim(sect_name)
        print "(a,i28)", ' # Cell section idx: ', sect
        print "(a,a28)", " # Bnd section type: ", &
          trim(ElementTypeName(cell_type))
        print "(a,i25)", " # Bnd condition color: ", &
          cgns_base(base) % block(block) % bnd_cond(bc) % color
        print "(a,i21)", " # Bnd section has # faces: ", &
          cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes
        print "(a,a25)", " # They belong to b.c.: ", &
          trim(cgns_base(base) % block(block) % bnd_cond(bc) % name)
        print "(a)", " #---------------------------------------------"
        print "(a)", " # Connections table (sample): "
        do loc = 1, min(6, cnt)
          print "(a,a16,4i8)", " # "," ", (face_n(n,loc), n = 1, n_nodes)
        end do
        print "(a)", " #---------------------------------------------"
        print "(a)", " # Parent data (sample): "
        !do loc = 1, min(6, cnt)
        !  print "(a,a16,4i8)", " # "," ", (parent_data(loc,n), n = 1, 2)
        !end do
      end if ! verbose

      do loc = 1, cnt
        ! Calculate absolute cell number
        cell = parent_data(loc, 1) + cnt_cells
        dir  = parent_data(loc, 2)
        grid % cells_bnd_color(dir, cell) = &
          cgns_base(base) % block(block) % bnd_cond(bc) % color

        if(loc < 7 .and. verbose) then
          print "(a,a24,3i8)", " # ", " ", cell, dir, &
            grid % cells_bnd_color(dir, cell)
        end if ! verbose

      end do ! loc = 1, cnt
    end if !bc_found .ne. 0
  end do ! bc = 1, cgns_base(base) % block(block) % n_bnd_conds

  deallocate(face_n)
  deallocate(parent_data)

  end subroutine
