!==============================================================================!
  subroutine Cgns_Mod_Read_2d_Bnd_Section_Connections_With_No_Parent_Data(  &
    base_id, block_id, sect_3d_id, grid, cell_n)
!------------------------------------------------------------------------------!
!   Read 2d elements connections info for current 3d sect and                  !
!   search for faces in given 3d cell connections list cell_n                  !
!   Apply b.c.                                                                 !
!   Prerequisites: all 3d cells in this block already follow                   !
!   grid structure https://cgns.github.io/CGNS_docs_current/sids/conv.html     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer              :: base_id, block_id, sect_3d_id
  integer, intent(in)  :: cell_n(:,:)
  type(Grid_Type)      :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base             ! base index number
  integer              :: block            ! block index number
  integer              :: sect_id          ! element section index
  character(len=80)    :: sect_name        ! name of the Elements_t node
  integer              :: cell_type        ! types of elements in the section
  integer              :: first_cell       ! index of first element
  integer              :: last_cell        ! index of last element
  integer, allocatable :: face_n(:,:)      ! 2d cells connections
  integer              :: n_nodes, loc, n, cell, dir, cnt, sect_3d
  integer              :: i, j, k, m, l, f_found, bc_found, sect, bc, f, m1
  logical              :: equal
  integer              :: error
  integer              :: parent_flag

  integer              :: cell_face(6,4)

  ! cgns HEXA_8 cell faces nodal connections
  integer, parameter, dimension(6, 4) :: face_hexa_8 =            &
                               transpose(reshape( (/ 1, 4, 3, 2,  &
                                                     1, 2, 6, 5,  &
                                                     2, 3, 7, 6,  &
                                                     3, 4, 8, 7,  &
                                                     1, 5, 8, 4,  &
                                                     5, 6, 7, 8  /), (/4, 6/) ))
  ! cgns PYRA_5 cell faces nodal connections
  integer, parameter, dimension(6, 4) :: face_pyra_5 =            &
                               transpose(reshape( (/ 1, 4, 3, 2,  &
                                                     1, 2, 5,-1,  &
                                                     2, 3, 5,-1,  &
                                                     3, 4, 5,-1,  &
                                                     4, 1, 5,-1,  &
                                                    -1,-1,-1,-1  /), (/4, 6/) ))
  ! cgns PENTA_6 cell faces nodal connections
  integer, parameter, dimension(6, 4) :: face_penta_6 =           &
                               transpose(reshape( (/ 1, 2, 5, 4,  &
                                                     2, 3, 6, 5,  &
                                                     3, 1, 4, 6,  &
                                                     1, 3, 2,-1,  &
                                                     4, 5, 6,-1,  &
                                                    -1,-1,-1,-1  /), (/4, 6/) ))
  ! cgns TETRA_4 cell faces nodal connections
  integer, parameter, dimension(6, 4) :: face_tetra_4 =           &
                               transpose(reshape( (/ 1, 3, 2,-1,  &
                                                     1, 2, 4,-1,  &
                                                     2, 3, 4,-1,  &
                                                     3, 1, 4,-1,  &
                                                    -1,-1,-1,-1,  &
                                                    -1,-1,-1,-1  /), (/4, 6/) ))
  ! 2d face
  integer, allocatable :: face_n_order(:,:)

  ! 3!
  integer, parameter, dimension(6, 3) :: face_3_n_order =      &
                               transpose(reshape( (/ 1, 2, 3,  &
                                                     1, 3, 2,  &
                                                     2, 1, 3,  &
                                                     2, 3, 1,  &
                                                     3, 1, 2,  &
                                                     3, 2, 1  /), (/3, 6/) ))

  ! 4!
  integer, parameter, dimension(24, 4) :: face_4_n_order =        &
                               transpose(reshape( (/ 1, 2, 3, 4,  &
                                                     2, 1, 3, 4,  &
                                                     3, 1, 2, 4,  &
                                                     1, 3, 2, 4,  &
                                                     2, 3, 1, 4,  &
                                                     3, 2, 1, 4,  &
                                                     3, 2, 4, 1,  &
                                                     2, 3, 4, 1,  &
                                                     4, 3, 2, 1,  &
                                                     3, 4, 2, 1,  &
                                                     2, 4, 3, 1,  &
                                                     4, 2, 3, 1,  &
                                                     4, 1, 3, 2,  &
                                                     1, 4, 3, 2,  &
                                                     3, 4, 1, 2,  &
                                                     4, 3, 1, 2,  &
                                                     1, 3, 4, 2,  &
                                                     3, 1, 4, 2,  &
                                                     2, 1, 4, 3,  &
                                                     1, 2, 4, 3,  &
                                                     4, 2, 1, 3,  &
                                                     2, 4, 1, 3,  &
                                                     1, 4, 2, 3,  &
                                                     4, 1, 2, 3 /), (/4, 24/) ))
!==============================================================================!

  ! Set input parameters
  base    = base_id
  block   = block_id
  sect_3d = sect_3d_id

  if ( size(cell_n, 1) .eq. 8 ) cell_face = face_hexa_8  ! HEXA_8
  if ( size(cell_n, 1) .eq. 5 ) cell_face = face_pyra_5  ! PYRA_5
  if ( size(cell_n, 1) .eq. 6 ) cell_face = face_penta_6 ! PENTA_6
  if ( size(cell_n, 1) .eq. 4 ) cell_face = face_tetra_4 ! TETRA_4

  do sect_id = 1, cgns_base(base) % block(block) % n_sects
    ! Introduce some abbreviations
    sect        = sect_id
    sect_name   = cgns_base(base) % block(block) % section(sect) % name
    cell_type   = cgns_base(base) % block(block) % section(sect) % cell_type
    first_cell  = cgns_base(base) % block(block) % section(sect) % first_cell
    last_cell   = cgns_base(base) % block(block) % section(sect) % last_cell
    parent_flag = cgns_base(base) % block(block) % section(sect) % parent_flag

    ! Number of cells in this section
    cnt = last_cell - first_cell + 1 ! cells in this sections

    ! Only 2d sections with no ParentData
    if ( ( ElementTypeName(cell_type) .eq. 'QUAD_4' .or.  &
           ElementTypeName(cell_type) .eq. 'TRI_3') .and. &
         (parent_flag .eq. 0)                          ) then

      ! Allocate memory
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') then
        n_nodes = 4
        allocate(face_n_order(24, n_nodes))
        face_n_order(:,:) = face_4_n_order(:,:)
      end if ! if ( ElementTypeName(cell_type) .eq. 'QUAD_4' ) then

      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) then
        n_nodes = 3
        allocate(face_n_order(6, n_nodes))
        face_n_order(:,:) = face_3_n_order(:,:)
      end if ! if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) then

      allocate(face_n(n_nodes, cnt))

      ! Read element data
      call Cg_Elements_Read_F(file_id,      & !(in )
                              base,         & !(in )
                              block,        & !(in )
                              sect,         & !(in )
                              face_n(:,:),  & !(out)
                              NULL,         & !(out) no parent data
                              error)          !(out)

      if (error.ne.0) then
        print "(a)", " #   Failed to read 2d elements connections from", sect_id
        call Cg_Error_Exit_F()
      endif

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
            print "(a,a24)", " #   2d cell section name: ", trim(sect_name)
            print "(a,a28)", " #   Bnd section type: ", &
              trim(ElementTypeName(cell_type))
            print "(a,i25)", " #   Bnd condition color: ", &
             cgns_base(base) % block(block) % bnd_cond(bc) % color
            print "(a,i21)", " #   Bnd section has # faces: ", &
              cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes
            print "(a,a25)", " #   They belong to b.c.: ", &
              trim(cgns_base(base) % block(block) % bnd_cond(bc) % name)
            print "(a)", " #-------------------------------------------------"
            print "(a)",     " #     Connections table (sample): "
            do loc = 1, min(6, cnt)
              print "(a,a16,4i8)", " # "," ", (face_n(n,loc), n = 1, n_nodes)
            end do
            print "(a)", " #-------------------------------------------------"
            print "(a)", " #   Recovered ParentData table (sample): "
          end if

          ! This loop should restore ParentData for current b.c. and sect_3d
          f_found = 0
          do j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes

            ! If point from b.c. list belongs to sect_3d
            equal = cgns_base(base) % block(block) % bnd_cond(bc) % &
              belongs_to_sect(j) .eq. sect_id
            if (.not. equal) cycle ! j loop

            ! Take face index corresponding to current bc from PointList
            f = cgns_base(base) % block(block) % bnd_cond(bc) % point_list(j)&
              - first_cell + 1 ! 1, cnt of face_n

            do i = 1, size(cell_n, 2)
              do k = 1, size(cell_face, 1) ! from 1 to 6th face
                do l = 1, size(face_n_order, 1) ! face nodal order (ex.,4-3-2-1)
                  m1 = 0
                  do m = 1, size(face_n_order, 2) ! from 1 to 3/4th node in face
                    if (cell_face(k,m).ne.-1) then
                      equal = &
                        face_n(face_n_order(l,m),f).eq.cell_n(cell_face(k,m),i)
                      if (.not. equal) exit ! m loop
                    else
                      m1 = m1 + 1
                      if ( m1 > 1 ) exit ! m loop
                    end if
                  end do ! m

                  ! If code reached this line, this means that 
                  ! 2d face connections are f_found inside 3d cell connections
                  ! with custom node order face_n_order(l,:) and this face
                  ! corresponds to bc since its index was taken from point_list

                  if (equal) then
                    f_found = f_found + 1
                    ! use Parent Data
                    ! like in Read_2d_Bnd_Section_Connections_With_Parent_Data
                    cell = cnt_cells + i
                    dir  = k
                    grid % cells_bnd_color(dir, cell) =  &
                         cgns_base(base) % block(block) % bnd_cond(bc) % color

                    if(verbose .and. f_found > 0 .and. f_found < 7) then
                      print "(a,a24,3i8)", " # ", " ", &
                        cell, dir, grid % cells_bnd_color(dir, cell)
                    end if

                  end if ! found 2d face in 3d cell
                end do ! l
              end do ! k
            end do ! i
          end do ! j

        if(verbose .and. f_found .ne. 0) then
          print "(a)",     " #-------------------------------------------------"
          print "(a,a22)", " #   Searching bnd faces in: ", &
            trim(cgns_base(base) % block(block) % section(sect_3d_id) % name)
          print "(a,i17)", " #   Found bnd faces in 3d block: ", f_found
          print "(a)",     " #-------------------------------------------------"
        end if ! verbose

      end if ! bc_f_found .ne. 0
    end do ! bc = 1, cgns_base(base) % block(block) % n_bnd_conds

    deallocate(face_n_order)
    deallocate(face_n)

    end if ! QUAD_4 and TRI_3
  end do ! 2d_sect = 1, cgns_base(base) % block(block) % n_sects

  end subroutine
