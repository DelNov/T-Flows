!==============================================================================!
  subroutine Load_Fluent(grid, file_name)
!------------------------------------------------------------------------------!
!   Reads the Fluent's file format.                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  character(SL)   :: file_name
!-----------------------------------[Locals]-----------------------------------!
  character(SL)        :: one_token
  integer              :: n_tri, n_quad, n_tet, n_hexa, n_pyra, n_wed, n_poly
  integer              :: n_cells, n_bnd_cells, n_faces, n_nodes
  integer              :: n_face_nodes, n_cells_zone
  integer              :: c, c1, c2, s, n, fu, i, l, pos
  integer              :: i_cel, i_nod, j_nod, k_nod, l_nod, i_pol
  integer              :: cell_type, face_type
  integer              :: cell_s, cell_e, side_s, side_e, node_s, node_e
  integer, allocatable :: cell_visited_from(:)
  integer              :: all_nodes(1024)      ! all cell's nodes
  integer              :: n_face_sect          ! number of face sections
  integer              :: face_sect_pos(2048)  ! where did Fluent store it
  integer              :: face_sect_bnd(2048)  ! where does T-Flows store it
  integer              :: n_bnd_cond           ! number of boundary conditions
  logical              :: this_sect_bnd        ! .true. if bnd cond section
  logical              :: the_end              ! end of file reached?
! character(SL)        :: bc_names(2048)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MIXED_ZONE = 0
  integer, parameter :: CELL_TRI   = 1
  integer, parameter :: CELL_TETRA = 2
  integer, parameter :: CELL_QUAD  = 3
  integer, parameter :: CELL_HEXA  = 4
  integer, parameter :: CELL_PYRA  = 5
  integer, parameter :: CELL_WEDGE = 6
  integer, parameter :: CELL_POLY  = 7
  integer, parameter :: FACE_TRI   = 3
  integer, parameter :: FACE_QUAD  = 4
!==============================================================================!

  call File_Mod_Open_File_For_Reading(file_name, fu)

  !-----------------------------------------------!
  !   Assume grid doesn't have polyhedral cells   !
  !-----------------------------------------------!
  grid % polyhedral = .false.

  !--------------------------------------------------------!
  !                                                        !
  !                                                        !
  !   Read number of nodes, cells, faces, boundary cells   !
  !                                                        !
  !                                                        !
  !--------------------------------------------------------!
  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Reading header                                           '
  print '(a60)', ' #----------------------------------------------------------'

  !------------------------------------------!
  !   Look for number of nodes and read it   !
  !     Nodes' section starts with '(10'     !
  !------------------------------------------!
  grid % n_nodes = 0
  rewind(fu)
  do while(grid % n_nodes .eq. 0)
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then
      if(line % tokens(1) .eq. '(10' .and. line % tokens(2) .eq. '(0') then
        read(line % tokens(4), '(z160)') grid % n_nodes
        print '(a34,i9)', ' # Number of nodes in header:     ', grid % n_nodes
      end if
    end if
  end do

  !------------------------------------------!
  !   Look for number of cells and read it   !
  !     Cells' section starts with '(12'     !
  !------------------------------------------!
  grid % n_cells = 0
  rewind(fu)
  do while(grid % n_cells .eq. 0)
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then
      if(line % tokens(1) .eq. '(12' .and. line % tokens(2) .eq. '(0') then
        read(line % tokens(4), '(z160)') grid % n_cells
        print '(a34,i9)', ' # Number of cells in header:     ', grid % n_cells
      end if
    end if
  end do

  !------------------------------------------!
  !   Look for number of faces and read it   !
  !     Faces' section starts with '(13'     !
  !------------------------------------------!
  grid % n_faces = 0
  rewind(fu)
  do while(grid % n_faces .eq. 0)
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then
      if(line % tokens(1) .eq. '(13' .and. line % tokens(2) .eq. '(0') then
        read(line % tokens(4), '(z160)') grid % n_faces
        print '(a34,i9)', ' # Number of faces in header:     ', grid % n_faces
      end if
    end if
  end do

  !----------------------------------------------------------------------------!
  !   Read the faces to count the boundary cells and boundary cells sections   !
  !----------------------------------------------------------------------------!
  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Reading face data to find number of boundary cells       '
  print '(a60)', ' #----------------------------------------------------------'
  n_face_sect = 0
  n_bnd_cond  = 0
  n_bnd_cells = 0
  n_faces     = 0
  rewind(fu)
  do while(n_faces < grid % n_faces)
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then

      ! Does the line mark the beginning of face-based data
      if(line % tokens(1) .eq. '(13' .and. line % tokens(2) .ne. '(0') then
        read(line % tokens(3), '(z160)') side_s  ! starting face
        read(line % tokens(4), '(z160)') side_e  ! ending face

        ! Increase the counter for face sections
        n_face_sect = n_face_sect + 1

        ! Store this face section position in the Fluent's mesh
        one_token = line % tokens(2)
        l = len_trim(one_token)
        read(one_token(2:l), '(z16)') face_sect_pos(n_face_sect)

        ! Take the cell type of this zone
        one_token = line % tokens(6)
        read(one_token(1:1), '(z1)') face_type
        if(face_type .eq. MIXED_ZONE) then
          print '(a34,i9,a4,i9)', ' # Found a mixed face zone from:  ',  &
                                  side_s, ' to:', side_e
        else
          print '(a34,i9,a4,i9)', ' # Found a uniform face zone from:',  &
                                  side_s, ' to:', side_e
        end if

        ! End the line if needed
        if(line % last .ne. '(') then
          call File_Mod_Read_Line(fu)
        end if

        !--------------------------!
        !   Browse through faces   !
        !--------------------------!
        this_sect_bnd = .false.
        do s = side_s, side_e
          n_faces = n_faces + 1
          call File_Mod_Read_Line(fu)

          ! Zone is mixed, read number of nodes and then c1 and c2
          if(face_type .eq. MIXED_ZONE) then
            read(line % tokens(1), *) n_face_nodes
            read(line % tokens(1+n_face_nodes+1), '(z160)') c1
            read(line % tokens(1+n_face_nodes+2), '(z160)') c2

          ! Zone is uniform, number of nodes is known, just read c1 and c2
          else
            if(face_type .eq. FACE_TRI)  n_face_nodes = 3
            if(face_type .eq. FACE_QUAD) n_face_nodes = 4
            read(line % tokens(0+n_face_nodes+1), '(z160)') c1
            read(line % tokens(0+n_face_nodes+2), '(z160)') c2
          end if

          ! This was a boundary face
          if( (c1 .eq. 0) .or. (c2 .eq. 0) ) then
            this_sect_bnd = .true.         ! marks this section as boundary
            n_bnd_cells = n_bnd_cells + 1
          end if
        end do

        !---------------------------------------------------------------!
        !   If this was a boundary section, store boundary conditions   !
        !---------------------------------------------------------------!
        if( this_sect_bnd ) then

          ! Increase the number of boundary condition sections and ...
          n_bnd_cond = n_bnd_cond + 1

          ! ... store mapping of this face section to T-Flows boundary condition
          face_sect_bnd(n_face_sect) = n_bnd_cond
        end if

      end if

    end if
  end do

  print '(a34,i9)', ' # Boundary cells from face data: ', n_bnd_cells
  print '(a34,i9)', ' # Boundary condition sections:   ', n_bnd_cond

  grid % n_bnd_cells = n_bnd_cells
  grid % n_bnd_cond  = n_bnd_cond
  allocate(grid % bnd_cond % name(n_bnd_cond))

  !--------------------------------------------!
  !                                            !
  !                                            !
  !   Allocate memory for Grid_Mod variables   !
  !                                            !
  !                                            !
  !--------------------------------------------!
  call Allocate_Memory(grid)

  !---------------------------!
  !                           !
  !                           !
  !   Read node coordinates   !
  !                           !
  !                           !
  !---------------------------!
  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Reading node coordinates                                 '
  print '(a60)', ' #----------------------------------------------------------'
  n_nodes = 0
  rewind(fu)
  do while(n_nodes < grid % n_nodes)
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then

      ! Does the line mark the beginning of node-based data
      if(line % tokens(1) .eq. '(10' .and. line % tokens(2) .ne. '(0') then
        read(line % tokens(3), '(z160)') node_s  ! starting node
        read(line % tokens(4), '(z160)') node_e  ! ending node
        print '(a34,i9,a4,i9)', ' # Found a node zone from:        ',  &
                                node_s, ' to:', node_e

        ! End the line if needed
        if(line % last .ne. '(') then
          call File_Mod_Read_Line(fu)
        end if

        ! Read all the nodes (node coordinates)
        do n = node_s, node_e
          n_nodes = n_nodes + 1
          call File_Mod_Read_Line(fu)
          read(line % tokens(1), *)  grid % xn(n)
          read(line % tokens(2), *)  grid % yn(n)
          read(line % tokens(3), *)  grid % zn(n)
        end do
      end if

    end if
  end do

  !-------------------------!
  !                         !
  !                         !
  !   Read the cell types   !
  !                         !
  !                         !
  !-------------------------!

  ! Initialize all cell counters
  n_cells = 0
  n_tri   = 0;  n_quad  = 0;  n_tet = 0
  n_hexa  = 0;  n_pyra  = 0;  n_wed = 0

  rewind(fu)
  do while(n_cells < grid % n_cells)
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then
      if(line % tokens(1) .eq. '(12' .and. line % tokens(2) .ne. '(0') then

        read(line % tokens(3), '(z160)') cell_s  ! starting cell
        read(line % tokens(4), '(z160)') cell_e  ! ending cell
        print '(a34,i9,a4,i9)', ' # Found a cell zone from:        ',  &
                                cell_s, ' to:', cell_e

        !---------------------------------------------------------!
        !   Check if the zone is mixed (listing all cell types)   !
        !---------------------------------------------------------!
        read(line % tokens(line % n_tokens)(1:1), '(z1)') cell_type

        !------------------------------------!
        !   You are reading a uniform zone   !
        !------------------------------------!
        if(cell_type .ne. MIXED_ZONE) then

          ! Number of cells in this zone
          n_cells_zone = cell_e - cell_s + 1

          ! Update the number of cells
          n_cells = n_cells + n_cells_zone

          if(cell_type .eq. CELL_TRI)   then
            n_tri = n_tri + n_cells_zone
            grid % cells_n_nodes(cell_s:cell_e) = 3
            grid % cells_n_faces(cell_s:cell_e) = 1
          end if
          if(cell_type .eq. CELL_QUAD)  then
            n_quad = n_quad + n_cells_zone
            grid % cells_n_nodes(cell_s:cell_e) = 4
            grid % cells_n_faces(cell_s:cell_e) = 1
          end if
          if(cell_type .eq. CELL_TETRA) then
            n_tet = n_tet + n_cells_zone
            grid % cells_n_nodes(cell_s:cell_e) = 4
            grid % cells_n_faces(cell_s:cell_e) = 4
          end if
          if(cell_type .eq. CELL_HEXA)  then
            n_hexa = n_hexa + n_cells_zone
            grid % cells_n_nodes(cell_s:cell_e) = 8
            grid % cells_n_faces(cell_s:cell_e) = 6
          end if
          if(cell_type .eq. CELL_PYRA)  then
            n_pyra  = n_pyra  + n_cells_zone
            grid % cells_n_nodes(cell_s:cell_e) = 5
            grid % cells_n_faces(cell_s:cell_e) = 5
          end if
          if(cell_type .eq. CELL_WEDGE) then
            n_wed = n_wed + n_cells_zone
            grid % cells_n_nodes(cell_s:cell_e) = 6
            grid % cells_n_faces(cell_s:cell_e) = 5
          end if

        !----------------------------------!
        !   You are reading a mixed zone   !
        !----------------------------------!
        else

6         continue

          call File_Mod_Read_Line(fu, remove='('//')')

          do i = 1, line % n_tokens
            read(line % tokens(i), *) cell_type

            ! Update the counter for all cells
            n_cells = n_cells + 1

            ! Update the counters for cell types
            if(cell_type .eq. CELL_TRI) then
              n_tri = n_tri + 1
              grid % cells_n_nodes(n_cells) = 3
              grid % cells_n_faces(n_cells) = 1

            else if(cell_type .eq. CELL_QUAD) then
              n_quad = n_quad + 1
              grid % cells_n_nodes(n_cells) = 4
              grid % cells_n_faces(n_cells) = 1

            else if(cell_type .eq. CELL_TETRA) then
              n_tet = n_tet + 1
              grid % cells_n_nodes(n_cells) = 4
              grid % cells_n_faces(n_cells) = 4

            else if(cell_type .eq. CELL_HEXA)  then
              n_hexa = n_hexa + 1
              grid % cells_n_nodes(n_cells) = 8
              grid % cells_n_faces(n_cells) = 6

            else if(cell_type .eq. CELL_PYRA)  then
              n_pyra = n_pyra + 1
              grid % cells_n_nodes(n_cells) = 5
              grid % cells_n_faces(n_cells) = 5

            else if(cell_type .eq. CELL_WEDGE) then
              n_wed = n_wed + 1
              grid % cells_n_nodes(n_cells) = 6
              grid % cells_n_faces(n_cells) = 5

            else if(cell_type .eq. CELL_POLY) then
              if(.not. grid % polyhedral) then
                print *, '# Found polyhedral cell(s), mesh is polyhedral!'
                grid % polyhedral = .true.
              end if
              n_poly = n_poly + 1
              grid % cells_n_nodes(n_cells) = -1
              grid % cells_n_faces(n_cells) = -1

            else
              print *, '# ERROR: Unsupported cell type', cell_type
              print *, '# This error is critical.  Exiting!'
              stop
            end if
          end do

          ! Did you reach the end of this secion?
          if(n_cells < cell_e) goto 6

        end if  ! end of the mixed zone, also end of both zones

      end if  ! found a beginning of cell zone
    end if  ! number of tokens bigger than one
  end do  ! infinite loop

  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Summary of interioir cell shapes:                        '
  print '(a60)', ' #----------------------------------------------------------'
  print '(a34,i9)', ' # Number of triangles:           ', n_tri
  print '(a34,i9)', ' # Number of quadrilaterals:      ', n_quad
  print '(a34,i9)', ' # Number of tetrahedra:          ', n_tet
  print '(a34,i9)', ' # Number of hexahedra:           ', n_hexa
  print '(a34,i9)', ' # Number of pyramids:            ', n_pyra
  print '(a34,i9)', ' # Number of wedges:              ', n_wed
  print '(a34,i9)', ' # Number of polyhedra:           ', n_poly

  !----------------------------------------------------!
  !                                                    !
  !                                                    !
  !   Read the faces again to store nodes, c1 and c2   !
  !                                                    !
  !                                                    !
  !----------------------------------------------------!
  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Reading face data to store c1, c2 and boundary cells     '
  print '(a60)', ' #----------------------------------------------------------'
  n_face_sect  = 0
  n_bnd_cells = 0
  n_faces     = 0
  rewind(fu)
  do while(n_faces < grid % n_faces)
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then

      ! Does the line mark the beginning of face-based data
      if(line % tokens(1) .eq. '(13' .and. line % tokens(2) .ne. '(0') then
        read(line % tokens(3), '(z160)') side_s  ! starting face
        read(line % tokens(4), '(z160)') side_e  ! ending face

        ! Increase the counter for face sections
        n_face_sect = n_face_sect + 1

        ! Take the cell type of this zone
        one_token = line % tokens(6)
        read(one_token(1:1), '(z1)') face_type
        if(face_type .eq. MIXED_ZONE) then
          print '(a34,i9,a4,i9)', ' # Found a mixed face zone from:  ',  &
                                  side_s, ' to:', side_e
        else
          print '(a34,i9,a4,i9)', ' # Found a uniform face zone from:',  &
                                  side_s, ' to:', side_e
        end if

        ! End the line if needed
        if(line % last .ne. '(') then
          call File_Mod_Read_Line(fu)
        end if

        !--------------------------!
        !   Browse through faces   !
        !--------------------------!
        do s = side_s, side_e
          n_faces = n_faces + 1
          call File_Mod_Read_Line(fu)

          ! Read and store number of nodes in this face
          if(face_type .eq. MIXED_ZONE) then
            read(line % tokens(1), *) n_face_nodes
          else
            if(face_type .eq. FACE_TRI)  n_face_nodes = 3
            if(face_type .eq. FACE_QUAD) n_face_nodes = 4
          end if
          grid % faces_n_nodes(s) = n_face_nodes

          ! Read nodes and cells surrounding the face for a mixed zone
          if(face_type .eq. MIXED_ZONE) then
            do i_nod = 1, n_face_nodes
              read(line % tokens(1+i_nod), '(z160)') grid % faces_n(i_nod, s)
            end do
            read(line % tokens(1+n_face_nodes+1), '(z160)') c1
            read(line % tokens(1+n_face_nodes+2), '(z160)') c2

          ! Read nodes and cells surrounding the face for a uniform zone
          else
            do i_nod = 1, n_face_nodes
              read(line % tokens(0+i_nod), '(z160)') grid % faces_n(i_nod, s)
            end do
            read(line % tokens(0+n_face_nodes+1), '(z160)') c1
            read(line % tokens(0+n_face_nodes+2), '(z160)') c2

          end if

          ! Case when c1 is a boundary cell
          if(c1 .eq. 0) then
            n_bnd_cells = n_bnd_cells + 1
            grid % faces_c(1, n_faces) = c2
            grid % faces_c(2, n_faces) = -n_bnd_cells
            grid % bnd_cond % color(-n_bnd_cells) = face_sect_bnd(n_face_sect)

          ! Case when c2 is a boundary cell
          else if(c2 .eq. 0) then
            n_bnd_cells = n_bnd_cells + 1
            grid % faces_c(1, n_faces) = c1
            grid % faces_c(2, n_faces) = -n_bnd_cells
            grid % bnd_cond % color(-n_bnd_cells) = face_sect_bnd(n_face_sect)

          ! Neither c1 nor c2 are boundary cells
          else
            grid % faces_c(1, n_faces) = min(c1, c2)
            grid % faces_c(2, n_faces) = max(c1, c2)

          end if

        end do
      end if

    end if
  end do

  print '(a34,i9)', ' # Boundary cells from face data: ', n_bnd_cells

  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Checking faces                                           '
  print '(a60)', ' #----------------------------------------------------------'

  !----------------------------------------!
  !   Check for duplicate nodes in faces   !
  !----------------------------------------!
  do s = 1, grid % n_faces
    n = grid % faces_n_nodes(s)
    do i_nod = 1, n
      do j_nod = i_nod+1, n
        if(grid % faces_n(i_nod, s) .eq. grid % faces_n(j_nod, s)) then
          print *, '# ERROR!  Duplicate nodes in face: ', s
          print *, '# This error is critical, exiting! '
          stop
        end if
      end do
    end do
  end do
  print *, '# No duplicate nodes in face data found, good!'

  !------------------------------!
  !                              !
  !                              !
  !   Reconstruct cells' nodes   !
  !                              !
  !                              !
  !------------------------------!
  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Reconstructing cells (determining their nodes)           '
  print '(a60)', ' #----------------------------------------------------------'
  allocate(cell_visited_from(grid % n_cells));  cell_visited_from(:) = 0

  !---------------------------------------------!
  !   Handle all boundary cells to start with   !
  !---------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    ! It is a boundary cell, just copy nodes
    if(c2 .lt. 0) then
      n = grid % faces_n_nodes(s)
      grid % cells_n_nodes(c2) = n
      grid % cells_n(1:n, c2)  = grid % faces_n(1:n, s)
    end if
  end do

  !-------------------------------------------------!
  !                                                 !
  !   First visits from faces; meaning that you:    !
  !                                                 !
  !   - browse through all quadrilateral faces to   !
  !     mark initial nodes for hexas and pyramids   !
  !                                                 !
  !   - browse through all triangular faces to      !
  !     mark initial nodes for tets and wedges      !
  !                                                 !
  !-------------------------------------------------!
  do s = 1, grid % n_faces

    !--------------------------------------!
    !   Handle cells c1 and c2 in one go   !
    !--------------------------------------!
    do i_cel = 1, 2

      c = grid % faces_c(i_cel, s)

      if(c .gt. 0) then

        !------------------------------------------------------------!
        !   First visit to pyramid or hexa from quadrilateral face   !
        !    For hexahedra that's face 5, for pyramid it's face 1    !
        !------------------------------------------------------------!
        if(grid % faces_n_nodes(s) .eq. 4) then
          if( (grid % cells_n_nodes(c) .eq. 5  .or.   &
               grid % cells_n_nodes(c) .eq. 8) .and.  &
              cell_visited_from(c) .eq. 0) then
            grid % cells_n(1:4, c) = grid % faces_n(1:4, s)
            cell_visited_from(c) = s
          end if
        end if

        !--------------------------------------------------------------!
        !   First visit to tetrahedron or wedge from triangular face   !
        !             For both shapes that would be face 1             !
        !--------------------------------------------------------------!
        if(grid % faces_n_nodes(s) .eq. 3) then
          if( (grid % cells_n_nodes(c) .eq. 4  .or.   &
               grid % cells_n_nodes(c) .eq. 6) .and.  &
              cell_visited_from(c) .eq. 0) then
            grid % cells_n(1:3, c) = grid % faces_n(1:3, s)
            cell_visited_from(c) = s
          end if
        end if  ! face is triangular

      end if  ! c > 0

    end do  ! through c1 and c2
  end do  ! through faces

  !---------------------------------------------------!
  !                                                   !
  !   Browse through all faces for the second visit   !
  !                                                   !
  !---------------------------------------------------!
  do s = 1, grid % n_faces

    !--------------------------------------!
    !   Handle cells c1 and c2 in one go   !
    !--------------------------------------!
    do i_cel = 1, 2

      c = grid % faces_c(i_cel, s)

      if(c .gt. 0) then

        !---------------------------!
        !   Face is quadrilateral   !
        !---------------------------!
        if(grid % faces_n_nodes(s) .eq. 4) then

          !--------------------------------!
          !   Second visit to hexahedron   !
          !--------------------------------!
          if( grid % cells_n_nodes(c) .eq. 8 .and.  &
              cell_visited_from(c) .ne.  s) then

            ! i_nod and j_nod, two consecutive nodes in the face
            do i_nod = 1, 4  ! nodes in face
              j_nod = i_nod + 1;  if(j_nod > 4) j_nod = j_nod - 4
              k_nod = i_nod + 2;  if(k_nod > 4) k_nod = k_nod - 4
              l_nod = i_nod + 3;  if(l_nod > 4) l_nod = l_nod - 4

              ! Face 1 same sense of rotation (see Cell_Numbering_Neu.f90)
              if(grid % faces_n(i_nod, s) .eq. grid % cells_n(1, c) .and.  &
                 grid % faces_n(j_nod, s) .eq. grid % cells_n(2, c) ) then
                grid % cells_n(6, c) = grid % faces_n(k_nod, s)
                grid % cells_n(5, c) = grid % faces_n(l_nod, s)
              end if

              ! Face 1 oposite sense of rotation (see Cell_Numbering_Neu.f90)
              if(grid % faces_n(i_nod, s) .eq. grid % cells_n(2, c) .and.  &
                 grid % faces_n(j_nod, s) .eq. grid % cells_n(1, c) ) then
                grid % cells_n(5, c) = grid % faces_n(k_nod, s)
                grid % cells_n(6, c) = grid % faces_n(l_nod, s)
              end if

              ! Face 3 same sense of rotation (see Cell_Numbering_Neu.f90)
              if(grid % faces_n(i_nod, s) .eq. grid % cells_n(3, c) .and.  &
                 grid % faces_n(j_nod, s) .eq. grid % cells_n(4, c) ) then
                grid % cells_n(8, c) = grid % faces_n(k_nod, s)
                grid % cells_n(7, c) = grid % faces_n(l_nod, s)
              end if

              ! Face 3 oposite sense of rotation (see Cell_Numbering_Neu.f90)
              if(grid % faces_n(i_nod, s) .eq. grid % cells_n(4, c) .and.  &
                 grid % faces_n(j_nod, s) .eq. grid % cells_n(3, c) ) then
                grid % cells_n(7, c) = grid % faces_n(k_nod, s)
                grid % cells_n(8, c) = grid % faces_n(l_nod, s)
              end if

            end do
          end if    ! if hexahedron

          !-----------------------------------------------------------!
          !   Second time you visit a wedge from quadrilateral face   !
          !-----------------------------------------------------------!
          if( grid % cells_n_nodes(c) .eq. 6 .and.  &
              cell_visited_from(c) .ne.  s) then

            ! i_nod and j_nod, two consecutive nodes in the face
            do i_nod = 1, 4  ! nodes in face
              j_nod = i_nod + 1;  if(j_nod > 4) j_nod = j_nod - 4
              k_nod = i_nod + 2;  if(k_nod > 4) k_nod = k_nod - 4
              l_nod = i_nod + 3;  if(l_nod > 4) l_nod = l_nod - 4

              ! Face 1 same sense of rotation (see Cell_Numbering_Neu.f90)
              if(grid % faces_n(i_nod, s) .eq. grid % cells_n(1, c) .and.  &
                 grid % faces_n(j_nod, s) .eq. grid % cells_n(2, c) ) then
                grid % cells_n(5, c) = grid % faces_n(k_nod, s)
                grid % cells_n(4, c) = grid % faces_n(l_nod, s)
              end if

              ! Face 1 oposite sense of rotation (see Cell_Numbering_Neu.f90)
              if(grid % faces_n(i_nod, s) .eq. grid % cells_n(2, c) .and.  &
                 grid % faces_n(j_nod, s) .eq. grid % cells_n(1, c) ) then
                grid % cells_n(4, c) = grid % faces_n(k_nod, s)
                grid % cells_n(5, c) = grid % faces_n(l_nod, s)
              end if

              ! Face 2 same sense of rotation (see Cell_Numbering_Neu.f90)
              if(grid % faces_n(i_nod, s) .eq. grid % cells_n(2, c) .and.  &
                 grid % faces_n(j_nod, s) .eq. grid % cells_n(3, c) ) then
                grid % cells_n(6, c) = grid % faces_n(k_nod, s)
                grid % cells_n(5, c) = grid % faces_n(l_nod, s)
              end if

              ! Face 2 oposite sense of rotation (see Cell_Numbering_Neu.f90)
              if(grid % faces_n(i_nod, s) .eq. grid % cells_n(3, c) .and.  &
                 grid % faces_n(j_nod, s) .eq. grid % cells_n(2, c) ) then
                grid % cells_n(5, c) = grid % faces_n(k_nod, s)
                grid % cells_n(6, c) = grid % faces_n(l_nod, s)
              end if

            end do

          end if

        end if  ! face is quadrilateral

        !------------------------!
        !   Face is triangular   !
        !------------------------!
        if(grid % faces_n_nodes(s) .eq. 3) then

          !-------------------------------------------------------------!
          !   Second time you visit a pyramid from triangular face      !
          !   (For pyramid, one node still missing, any face will do)   !
          !-------------------------------------------------------------!
          if( grid % cells_n_nodes(c) .eq. 5 .and.  &
              cell_visited_from(c) .ne.  s   .and.  &
              cell_visited_from(c) .ne. -1) then
            do i_nod = 1, 3
              if(all(grid % cells_n(1:4, c) .ne.  &
                     grid % faces_n(i_nod, s))) then
                grid % cells_n(5, c) = grid % faces_n(i_nod, s)
              end if
            end do

            ! Don't visit this cell again
            cell_visited_from(c) = -1
          end if

          !-----------------------------------------------------------------!
          !   Second time you visit a tetrahedron from triangular face      !
          !   (For tetrahedron, one node still missing, any face will do)   !
          !-----------------------------------------------------------------!
          if( grid % cells_n_nodes(c) .eq. 4 .and.  &
              cell_visited_from(c) .ne.  s   .and.  &
              cell_visited_from(c) .ne. -1) then
            do i_nod = 1, 3
              if(all(grid % cells_n(1:3, c) .ne.  &
                     grid % faces_n(i_nod, s))) then
                grid % cells_n(4, c) = grid % faces_n(i_nod, s)
              end if
            end do

            ! Don't visit this cell again
            cell_visited_from(c) = -1
          end if
        end if  ! face is triangular

      end if  ! c > 0

    end do  ! through c1 and c2
  end do  ! through s, faces

  !----------------------------------------!
  !   Check for duplicate nodes in cells   !
  !----------------------------------------!
  do c = 1, grid % n_cells
    n = grid % cells_n_nodes(c)
    do i_nod = 1, n
      do j_nod = i_nod+1, n
        if(grid % cells_n(i_nod, c) .eq. grid % cells_n(j_nod, c)) then
          print *, '# ERROR!  Duplicate nodes in cell: ', c
          print *, '# This error is critical, exiting! '
          stop
        end if
      end do
    end do
  end do
  print *, '# No duplicate nodes in cell data found, good!'

  !------------------------------------------------------------!
  !                                                            !
  !   Count the faces of all polyhedral cells and store them   !
  !                                                            !
  !------------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    if( grid % cells_n_nodes(c1) .eq. -1) then
      if( cell_visited_from(c1) .ne. 0 ) then
        print *, '# A strange error has occurred!'
      end if

      ! Increase the number of polygonal faces for c1
      grid % cells_n_polyf(c1) = grid % cells_n_polyf(c1) + 1
      n = grid % cells_n_polyf(c1)
      grid % cells_p(n, c1) = s
    end if

    if( grid % cells_n_nodes(c2) .eq. -1) then
      if( cell_visited_from(c2) .ne. 0 ) then
        print *, '# A strange error has occurred!'
      end if

      ! Increase the number of polygonal faces for c2 and store the face
      grid % cells_n_polyf(c2) = grid % cells_n_polyf(c2) + 1
      n = grid % cells_n_polyf(c2)
      grid % cells_p(n, c2) = s
    end if

  end do

  !----------------------------------------------------------------!
  !                                                                !
  !   With faces counted, you can also store nodes for each cell   !
  !                                                                !
  !----------------------------------------------------------------!
  do c = 1, grid % n_cells

    ! Only do this for polyhedral cells
    if(grid % cells_n_nodes(c) .eq. -1) then

      ! Accumulate nodes from all faces surrounding the cell
      n = 0
      do i_pol = 1, grid % cells_n_polyf(c)
        s = grid % cells_p(i_pol, c)           ! take true face index
        do i_nod = 1, grid % faces_n_nodes(s)
          n = n + 1
          all_nodes(n) = grid % faces_n(i_nod, s)
        end do
      end do

      ! Perform a unique sort to remove duplicates
      call Sort_Mod_Unique_Int(all_nodes(1:n), n)

      grid % cells_n(1:n, c)  = all_nodes(1:n)
      grid % cells_n_nodes(c) = -n

    end if  ! if cell was polyhedral

  end do

  deallocate(cell_visited_from)

  !------------------------------------------------------------!
  !                                                            !
  !   Browse through all sections (fluid, boundary, and what   !
  !    ever) to store number of boundary conditions to grid    !
  !                                                            !
  !------------------------------------------------------------!

  DO N = 1, 8
    PRINT *, 'FACE_SECT_POS(N), FACE_SECT_BND(N)',  &
              FACE_SECT_POS(N), FACE_SECT_BND(N)
  END DO

  rewind(fu)
  the_end = .false.
  do while(.not. the_end)
    call File_Mod_Read_Line(fu, reached_end=the_end)
    if(.not. the_end .and. line % n_tokens > 1) then
      if(line % tokens(1) .eq. '(39' .or.  &
         line % tokens(1) .eq. '(45') then

        PRINT *, 'FOUND A SECTION!'
        PRINT *, 'WHOLE LINE: ', LINE % WHOLE

        ! Extract this face section position in the Fluent's mesh
        one_token = line % tokens(2)
        l = len_trim(one_token)
        read(one_token(2:l), '(i16)') pos
        PRINT *, '  THIS SECTION''S POSITION: ', POS

        ! Extract boundary condition name
        one_token = line % tokens(4)(1:index(line % tokens(4), ')')-1)

        do n = 1, 2048
          if(face_sect_bnd(n) .ne. 0) then
            if(face_sect_pos(n) .eq. pos) then
              PRINT *, '    THIS SECTION IS A BOUNDARY SECTION'
              PRINT *, '    THIS SECTION BC NUMBER IS: ', FACE_SECT_BND(N)
              grid % bnd_cond % name(face_sect_bnd(n)) = one_token
              PRINT *, '    THIS SECTION BC NAME IS    ',  &
                       GRID % BND_COND % NAME(FACE_SECT_BND(N))
            end if
          end if
        end do

      end if
    end if
  end do
  STOP

  close(fu)

  end subroutine
