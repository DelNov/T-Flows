!==============================================================================!
  subroutine Load_Fluent(Convert, Grid, file_name)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to read grid files in the Fluent's file format
!>  and populate the provided Grid object with the necessary data.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * File opening and initial setup: The subroutine opens the Fluent file in  !
!     binary mode, assuming that the grid does not contain polyhedral cells    !
!     initially.                                                               !
!   * Reading header information: It extracts essential data from the file     !
!     header, including the number of nodes, cells, and faces in the mesh.     !
!   * Processing face data: The subroutine counts boundary cells and detrmines !
!     the number of boundary cell sections by reading through the face data.   !
!     This includes identifying the type of each face (mixed, quadrilateral,   !
!     or triangular) and the cells (c1 and c2) adjacent to each face.          !
!   * Memory allocation: Based on the information gathered, it allocates       !
!     memory  for various Grid structure variables to store the mesh data.     !
!   * Reading node coordinates: The subroutine reads and stores the            !
!     coordinates of each node in the grid.                                    !
!   * Reading and classifying cells: It reads through cell data to classify    !
!     cells into different types (e.g., triangles, quadrilaterals, tetrahedra, !
!     hexahedra, pyramids, wedges, polyhedra) and counts the number of each    !
!     type.                                                                    !
!   * Processing faces again for detailed data: The subroutine revisits the    !
!     face data to store detailed information about each face, including nodes !
!     and adjacent cells.
!   * Reconstructing cell nodes: For non-polyhedral cells, it reconstructs the !
!     nodes of each cell based on face data. For polyhedral cells, nodes are   !
!     accumulated from all surrounding faces.                                  !
!   * Data verification and integrity checks: It performs various checks to    !
!     ensure the integrity of the data, such as verifying that no duplicate    !
!     nodes or cells exist in the face and cell data.                          !
!   * Final Processing of Boundary Conditions: The subroutine goes through all !
!     sections in the file to store the number of boundary conditions and      !
!     their names in the Grid structure.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert    !! parent class
  type(Grid_Type)     :: Grid       !! grid being converted
  character(SL)       :: file_name  !! file name
!----------------------------------[Calling]-----------------------------------!
# ifdef __NVCOMPILER_LLVM__
  integer :: ftell  ! needed for Nvidia compiler
# endif
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
  integer, parameter :: FACE_POLY  = 5  ! an educated guess, but seems good
!-----------------------------------[Locals]-----------------------------------!
  character(SL)          :: one_token
  character(1)           :: one_char
  integer(DP)            :: offset
  integer                :: n_tri, n_quad, n_tet, n_hexa, n_pyra, n_wed, n_poly
  integer                :: n_cells, n_bnd_cells, n_faces, n_nodes
  integer                :: n_face_nodes, n_cells_zone
  integer                :: c, c1, c2, s, n, fu, l, pos, length, error
  integer                :: i_cel, i_nod, j_nod, k_nod, l_nod, i_fac
  integer                :: cell_type, zone_type
  integer                :: cell_f, cell_l, side_f, side_l, node_f, node_l
  integer                :: all_nodes(1024)       ! all cell's nodes
  integer                :: n_face_sect           ! number of face sections
  integer                :: face_sect_pos(2048)   ! where did Fluent store it
  integer                :: face_sect_bnd(2048)   ! where does T-Flows store it
  integer                :: n_bnd_reg            ! number of boundary conditions
  logical                :: this_sect_bnd         ! .true. if bnd cond section
  logical                :: the_end               ! end of file reached?
  logical                :: ascii                 ! is file in ascii format?
  integer,   allocatable :: cell_visited_from(:), cell_types(:)
!==============================================================================!

  call Profiler % Start('Load_Fluent')

  ! Open the file in binary mode, because it could be mixed
  call File % Open_For_Reading_Binary(file_name, fu)

  !-----------------------------------------------!
  !   Assume grid doesn't have polyhedral cells   !
  !-----------------------------------------------!
  Grid % polyhedral = .false.

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
  !                                          !
  !   Look for number of nodes and read it   !
  !     Nodes' section starts with '(10'     !
  !                                          !
  !------------------------------------------!
  Grid % n_nodes = 0
  rewind(fu)
  do while(Grid % n_nodes .eq. 0)
    call File % Read_Line(fu)
    if(Line % n_tokens > 1) then
      if(Line % tokens(1) .eq. '(10' .and. Line % tokens(2) .eq. '(0') then
        read(Line % tokens(4), '(z160)') Grid % n_nodes
        print '(a34,i9)', ' # Number of nodes in header:     ', Grid % n_nodes
      end if
    end if
  end do

  !------------------------------------------!
  !                                          !
  !   Look for number of cells and read it   !
  !     Cells' section starts with '(12'     !
  !                                          !
  !------------------------------------------!
  Grid % n_cells = 0
  rewind(fu)
  do while(Grid % n_cells .eq. 0)
    call File % Read_Line(fu)
    if(Line % n_tokens > 1) then
      if(Line % tokens(1) .eq. '(12' .and. Line % tokens(2) .eq. '(0') then
        read(Line % tokens(4), '(z160)') Grid % n_cells
        print '(a34,i9)', ' # Number of cells in header:     ', Grid % n_cells
      end if
    end if
  end do

  !------------------------------------------!
  !                                          !
  !   Look for number of faces and read it   !
  !     Faces' section starts with '(13'     !
  !                                          !
  !------------------------------------------!
  Grid % n_faces = 0
  rewind(fu)
  do while(Grid % n_faces .eq. 0)
    call File % Read_Line(fu)
    if(Line % n_tokens > 1) then
      if(Line % tokens(1) .eq. '(13' .and. Line % tokens(2) .eq. '(0') then
        read(Line % tokens(4), '(z160)') Grid % n_faces
        print '(a34,i9)', ' # Number of faces in header:     ', Grid % n_faces
      end if
    end if
  end do

  !----------------------------------------------------------------------------!
  !                                                                            !
  !   Read the faces to count the boundary cells and boundary cells sections   !
  !                                                                            !
  !----------------------------------------------------------------------------!
  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Reading face data to find number of boundary cells       '
  print '(a60)', ' #----------------------------------------------------------'
  face_sect_pos(:) = 0
  face_sect_bnd(:) = 0
  n_face_sect      = 0
  n_bnd_reg        = 0
  n_bnd_cells      = 0
  n_faces          = 0
  the_end          = .false.
  rewind(fu)
  do while(n_faces < Grid % n_faces .and. .not. the_end)
    call File % Read_Line(fu, reached_end=the_end)
    if(Line % n_tokens > 1) then

      !----------------------------------------------------------!
      !   Does the line mark the beginning of face-based data?   !
      !   13) is for ascii and 2013) for binary formatted data   !
      !----------------------------------------------------------!
      if( (Line % tokens(1) .eq. '(13' .or. Line % tokens(1) .eq. '(2013')  &
         .and. Line % tokens(2) .ne. '(0') then

        ! Store the format of this section
        ascii = .true.
        if(Line % tokens(1) .eq. '(2013') ascii = .false.

        ! Fetch start and ending face
        read(Line % tokens(3), '(z160)') side_f  ! first face
        read(Line % tokens(4), '(z160)') side_l  ! last face

        ! Increase the counter for face sections
        n_face_sect = n_face_sect + 1

        ! Store this face section position in the Fluent's mesh
        one_token = Line % tokens(2)
        l = len_trim(one_token)
        read(one_token(2:l), '(z16)') face_sect_pos(n_face_sect)

        ! Take the cell type of this zone
        one_token = Line % tokens(6)
        read(one_token(1:1), '(z1)') zone_type
        if(zone_type .eq. MIXED_ZONE) then
          print '(a34,i9,a4,i9)', ' # Found a mixed face zone from:  ',  &
                                  side_f, ' to:', side_l
        else
          print '(a34,i9,a4,i9)', ' # Found a uniform face zone from:',  &
                                  side_f, ' to:', side_l
        end if

        ! End the line if needed, just read one left bracket '('
        if(Line % last .ne. '(') then
          if(ascii)       call File % Read_Line(fu)
          if(.not. ascii) read(fu) one_char
        end if

        !--------------------------!
        !   Browse through faces   !
        !--------------------------!
        this_sect_bnd = .false.
        do s = side_f, side_l
          n_faces = n_faces + 1

          ! Zone is mixed, read number of nodes and then c1 and c2
          if(zone_type .eq. MIXED_ZONE .or.  &
             zone_type .eq. FACE_POLY) then
            if(ascii) then
              call File % Read_Line(fu)
              read(Line % tokens(1), *) n_face_nodes
              read(Line % tokens(1+n_face_nodes+1), '(z160)') c1
              read(Line % tokens(1+n_face_nodes+2), '(z160)') c2
            else
              call File % Read_Binary_Int4_Array(fu, 1)
              n_face_nodes = int4_array(1)
              call File % Read_Binary_Int4_Array(fu, n_face_nodes)
              call File % Read_Binary_Int4_Array(fu, 2)
              c1 = int4_array(1)
              c2 = int4_array(2)
            end if

          ! Zone is uniform, number of nodes is known, just read c1 and c2
          else
            if(zone_type .eq. FACE_TRI)  n_face_nodes = 3
            if(zone_type .eq. FACE_QUAD) n_face_nodes = 4
            if(ascii) then
              call File % Read_Line(fu)
              read(Line % tokens(0+n_face_nodes+1), '(z160)') c1
              read(Line % tokens(0+n_face_nodes+2), '(z160)') c2
            else
              call File % Read_Binary_Int4_Array(fu, n_face_nodes)
              call File % Read_Binary_Int4_Array(fu, 2)
              c1 = int4_array(1)
              c2 = int4_array(2)
            end if
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
          n_bnd_reg = n_bnd_reg + 1

          ! ... store mapping of this face section to T-Flows boundary condition
          face_sect_bnd(n_face_sect) = n_bnd_reg
        end if

      end if

    end if
  end do
  if(the_end) then
    print *, '# ERROR: Could not find face-based data in: ', file_name
    print *, '# This error is critical!  Exiting!'
    stop
  end if

  print '(a34,i9)', ' # Boundary cells from face data: ', n_bnd_cells
  print '(a34,i9)', ' # Boundary condition sections:   ', n_bnd_reg

  Grid % n_bnd_cells = n_bnd_cells
  call Grid % Allocate_Regions(n_bnd_reg)

  !--------------------------------------------!
  !                                            !
  !                                            !
  !   Allocate memory for Grid_Mod variables   !
  !                                            !
  !                                            !
  !--------------------------------------------!
  call Convert % Allocate_Memory(Grid)

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
  the_end = .false.
  rewind(fu)
  do while(n_nodes < Grid % n_nodes .and. .not. the_end)
    call File % Read_Line(fu)
    if(Line % n_tokens > 1) then

      !----------------------------------------------------------!
      !   Does the line mark the beginning of node-based data?   !
      !   10) is for ascii and 3010) for binary formatted data   !
      !----------------------------------------------------------!
      if( (Line % tokens(1) .eq. '(10' .or. Line % tokens(1) .eq. '(3010')  &
         .and. Line % tokens(2) .ne. '(0') then
        read(Line % tokens(3), '(z160)') node_f  ! first node
        read(Line % tokens(4), '(z160)') node_l  ! last node
        print '(a34,i9,a4,i9)', ' # Found a node zone from:        ',  &
                                node_f, ' to:', node_l

        ! Store the format of this section
        ascii = .true.
        if(Line % tokens(1) .eq. '(3010') ascii = .false.

        ! End the Line if needed
        if(Line % last .ne. '(') then
          if(ascii)       call File % Read_Line(fu)
          if(.not. ascii) read(fu) one_char
        end if

        ! Read all the nodes (node coordinates)
        do n = node_f, node_l
          n_nodes = n_nodes + 1
          if(ascii) then
            call File % Read_Line(fu)
            read(Line % tokens(1), *)  Grid % xn(n)
            read(Line % tokens(2), *)  Grid % yn(n)
            read(Line % tokens(3), *)  Grid % zn(n)
          else
            call File % Read_Binary_Real8_Array(fu, 3)
            Grid % xn(n) = real8_array(1)
            Grid % yn(n) = real8_array(2)
            Grid % zn(n) = real8_array(3)
          end if
        end do

        if(.not. ascii) read(fu) one_char
      end if

    end if
  end do
  if(the_end) then
    print *, '# ERROR: Could not find node-based data in: ', file_name
    print *, '# This error is critical!  Exiting!'
    stop
  end if

  !-------------------------!
  !                         !
  !                         !
  !   Read the cell types   !
  !                         !
  !                         !
  !-------------------------!

  ! Initialize all cell counters
  n_cells = 0
  n_tri   = 0
  n_quad  = 0
  n_tet   = 0
  n_hexa  = 0
  n_pyra  = 0
  n_wed   = 0
  n_poly  = 0
  the_end = .false.

  rewind(fu)
  do while(n_cells < Grid % n_cells .and. .not. the_end)
    call File % Read_Line(fu, reached_end=the_end)
    if(Line % n_tokens > 1) then
      if( (Line % tokens(1) .eq. '(12' .or. Line % tokens(1) .eq. '(2012')  &
         .and. Line % tokens(2) .ne. '(0') then

        ! Store the format of this section
        ascii = .true.
        if(Line % tokens(1) .eq. '(2012') ascii = .false.

        ! Fetch first and last cell
        read(Line % tokens(3), '(z160)') cell_f  ! first cell
        read(Line % tokens(4), '(z160)') cell_l  ! last cell

        ! Check if the zone is mixed (listing all cell types)
        read(Line % tokens(6)(1:1), '(z1)') zone_type
        if(zone_type .eq. MIXED_ZONE) then
          print '(a34,i9,a4,i9)', ' # Found a mixed cell zone from:  ',  &
                                  cell_f, ' to:', cell_l
        else
          print '(a34,i9,a4,i9)', ' # Found a uniform cell zone from:',  &
                                  cell_f, ' to:', cell_l
        end if
        if(.not. ascii) read(fu) one_char  ! read "("

        !------------------------------------!
        !   You are reading a uniform zone   !
        !------------------------------------!
        if(zone_type .ne. MIXED_ZONE) then

          ! Number of cells in this zone
          n_cells_zone = cell_l - cell_f + 1

          ! Update the number of cells
          n_cells = n_cells + n_cells_zone

          if(zone_type .eq. CELL_TRI)   then
            n_tri = n_tri + n_cells_zone
            Grid % cells_n_nodes(cell_f:cell_l) = 3
            Grid % cells_n_faces(cell_f:cell_l) = 1
          end if
          if(zone_type .eq. CELL_QUAD)  then
            n_quad = n_quad + n_cells_zone
            Grid % cells_n_nodes(cell_f:cell_l) = 4
            Grid % cells_n_faces(cell_f:cell_l) = 1
          end if
          if(zone_type .eq. CELL_TETRA) then
            n_tet = n_tet + n_cells_zone
            Grid % cells_n_nodes(cell_f:cell_l) = 4
            Grid % cells_n_faces(cell_f:cell_l) = 4
          end if
          if(zone_type .eq. CELL_HEXA)  then
            n_hexa = n_hexa + n_cells_zone
            Grid % cells_n_nodes(cell_f:cell_l) = 8
            Grid % cells_n_faces(cell_f:cell_l) = 6
          end if
          if(zone_type .eq. CELL_PYRA)  then
            n_pyra  = n_pyra  + n_cells_zone
            Grid % cells_n_nodes(cell_f:cell_l) = 5
            Grid % cells_n_faces(cell_f:cell_l) = 5
          end if
          if(zone_type .eq. CELL_WEDGE) then
            n_wed = n_wed + n_cells_zone
            Grid % cells_n_nodes(cell_f:cell_l) = 6
            Grid % cells_n_faces(cell_f:cell_l) = 5
          end if
          if(zone_type .eq. CELL_POLY)  then
            n_poly  = n_poly  + n_cells_zone
            Grid % cells_n_nodes(cell_f:cell_l) = -1  ! attend
            Grid % cells_n_faces(cell_f:cell_l) = -1  ! have to attend too
            Grid % polyhedral = .true.
          end if

        !----------------------------------!
        !   You are reading a mixed zone   !
        !----------------------------------!
        else

          n_cells_zone = cell_l - cell_f + 1  ! this was read above

          if(ascii) then
            ! Find out the cell record length
            ! This is approximate, just until you encounter next
            ! ')' character, but will nonetheless do the job
            offset = ftell(fu)  ! mark the offset
            length = 0
            do
              read(fu) one_char
              if(one_char .eq. ')') exit
              length = length + 1
            end do

            ! Once you've found lenght, go back to the beginning
            ! of the record using the offset you stored before
#           ifdef __INTEL_COMPILER
            error = fseek(fu, offset, 0)
#           else
            error = 0
            call fseek(fu, offset, 0)
#           endif

            ! Allocate helping arrays
            allocate(cell_types(n_cells_zone))  ! allocate cell types

            ! Browse the very long Line and read cell types from it
            c = 0  ! cell counter
            do
              read(fu) one_char
              if(one_char .eq. ')') exit
              if(one_char .ge. '0' .and. one_char .le. '9') then
                c = c + 1                             ! increase cell count
                read(one_char, '(i1)') cell_types(c)  ! read cell types
              end if
            end do

            ! At this point: c .eq. n_cells_zone: check it!
            Assert(c .eq. n_cells_zone)
          end if

          do c = 1, n_cells_zone

            if(ascii) then
              cell_type = cell_types(c)
            else
              call File % Read_Binary_Int4_Array(fu, 1)
              cell_type = int4_array(1)
            end if

            ! Update the counter for all cells
            n_cells = n_cells + 1

            ! Update the counters for cell types
            if(cell_type .eq. CELL_TRI) then
              n_tri = n_tri + 1
              Grid % cells_n_nodes(n_cells) = 3
              Grid % cells_n_faces(n_cells) = 1

            else if(cell_type .eq. CELL_QUAD) then
              n_quad = n_quad + 1
              Grid % cells_n_nodes(n_cells) = 4
              Grid % cells_n_faces(n_cells) = 1

            else if(cell_type .eq. CELL_TETRA) then
              n_tet = n_tet + 1
              Grid % cells_n_nodes(n_cells) = 4
              Grid % cells_n_faces(n_cells) = 4

            else if(cell_type .eq. CELL_HEXA)  then
              n_hexa = n_hexa + 1
              Grid % cells_n_nodes(n_cells) = 8
              Grid % cells_n_faces(n_cells) = 6

            else if(cell_type .eq. CELL_PYRA)  then
              n_pyra = n_pyra + 1
              Grid % cells_n_nodes(n_cells) = 5
              Grid % cells_n_faces(n_cells) = 5

            else if(cell_type .eq. CELL_WEDGE) then
              n_wed = n_wed + 1
              Grid % cells_n_nodes(n_cells) = 6
              Grid % cells_n_faces(n_cells) = 5

            else if(cell_type .eq. CELL_POLY) then
              if(.not. Grid % polyhedral) then
                print *, '# Found polyhedral cell(s), mesh is polyhedral!'
                Grid % polyhedral = .true.
              end if
              n_poly = n_poly + 1
              Grid % cells_n_nodes(n_cells) = -1
              Grid % cells_n_faces(n_cells) = -1

            else
              print *, '# ERROR: Unsupported cell type', cell_type
              print *, '# This error is critical.  Exiting!'
              stop
            end if
          end do

        end if  ! end of the mixed zone, also end of both zones

      end if  ! found a beginning of cell zone
    end if  ! number of tokens bigger than one
  end do  ! infinite loop
  if(the_end) then
    print *, '# ERROR: Could not find cell-based data in: ', file_name
    print *, '# This error is critical!  Exiting!'
    stop
  end if

  ! This looks like a convenient time to adjust the dimenson for cells' nodes
  call Adjust_First_Dim(maxval(Grid % cells_n_nodes(1:n_cells)), Grid % cells_n)

  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Summary of interior cell shapes:                         '
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
  the_end     = .false.
  rewind(fu)
  do while(n_faces < Grid % n_faces .and. .not. the_end)
    call File % Read_Line(fu, reached_end=the_end)
    if(Line % n_tokens > 1) then

      !----------------------------------------------------------!
      !   Does the Line mark the beginning of face-based data?   !
      !   13) is for ascii and 2013) for binary formatted data   !
      !----------------------------------------------------------!
      if( (Line % tokens(1) .eq. '(13' .or. Line % tokens(1) .eq. '(2013')  &
         .and. Line % tokens(2) .ne. '(0') then

        ! Store the format of this section
        ascii = .true.
        if(Line % tokens(1) .eq. '(2013') ascii = .false.

        ! Fetch first and last face
        read(Line % tokens(3), '(z160)') side_f  ! first face
        read(Line % tokens(4), '(z160)') side_l  ! last face

        ! Increase the counter for face sections
        n_face_sect = n_face_sect + 1

        ! Take the cell type of this zone
        one_token = Line % tokens(6)
        read(one_token(1:1), '(z1)') zone_type
        if(zone_type .eq. MIXED_ZONE) then
          print '(a34,i9,a4,i9)', ' # Found a mixed face zone from:  ',  &
                                  side_f, ' to:', side_l
        else
          print '(a34,i9,a4,i9)', ' # Found a uniform face zone from:',  &
                                  side_f, ' to:', side_l
        end if

        ! End the Line if needed
        if(Line % last .ne. '(') then
          if(ascii)       call File % Read_Line(fu)
          if(.not. ascii) read(fu) one_char
        end if

        !--------------------------!
        !   Browse through faces   !
        !--------------------------!
        do s = side_f, side_l
          n_faces = n_faces + 1

          ! Read nodes and cells surrounding the face for a mixed zone
          if(zone_type .eq. MIXED_ZONE .or.  &
             zone_type .eq. FACE_POLY) then
            if(ascii) then
              call File % Read_Line(fu)
              read(Line % tokens(1), *) n_face_nodes
              Grid % faces_n_nodes(s) = n_face_nodes
              call Adjust_First_Dim(n_face_nodes, Grid % faces_n)
              do i_nod = 1, n_face_nodes
                read(Line % tokens(1+i_nod), '(z160)') Grid % faces_n(i_nod, s)
              end do
              read(Line % tokens(1+n_face_nodes+1), '(z160)') c1
              read(Line % tokens(1+n_face_nodes+2), '(z160)') c2
            else
              call File % Read_Binary_Int4_Array(fu, 1)
              n_face_nodes = int4_array(1)
              call File % Read_Binary_Int4_Array(fu, n_face_nodes)
              call Adjust_First_Dim(n_face_nodes, Grid % faces_n)
              do i_nod = 1, n_face_nodes
                Grid % faces_n(i_nod, s) = int4_array(i_nod)
              end do
              call File % Read_Binary_Int4_Array(fu, 2)
              c1 = int4_array(1)
              c2 = int4_array(2)
            end if

          ! Read nodes and cells surrounding the face for a uniform zone
          else
            if(zone_type .eq. FACE_TRI)  n_face_nodes = 3
            if(zone_type .eq. FACE_QUAD) n_face_nodes = 4
            call Adjust_First_Dim(n_face_nodes, Grid % faces_n)
            if(ascii) then
              call File % Read_Line(fu)
              do i_nod = 1, n_face_nodes
                read(Line % tokens(0+i_nod), '(z160)') Grid % faces_n(i_nod, s)
              end do
              read(Line % tokens(0+n_face_nodes+1), '(z160)') c1
              read(Line % tokens(0+n_face_nodes+2), '(z160)') c2
            else
              call File % Read_Binary_Int4_Array(fu, n_face_nodes)
              do i_nod = 1, n_face_nodes
                Grid % faces_n(i_nod, s) = int4_array(i_nod)
              end do
              call File % Read_Binary_Int4_Array(fu, 2)
              c1 = int4_array(1)
              c2 = int4_array(2)
            end if

          end if
          Grid % faces_n_nodes(s) = n_face_nodes

          ! Case when c1 is a boundary cell
          if(c1 .eq. 0) then
            n_bnd_cells = n_bnd_cells + 1
            Grid % faces_c(1, s) = c2
            Grid % faces_c(2, s) = -n_bnd_cells
            Grid % region % at_cell(-n_bnd_cells) = face_sect_bnd(n_face_sect)

          ! Case when c2 is a boundary cell
          else if(c2 .eq. 0) then
            n_bnd_cells = n_bnd_cells + 1
            Grid % faces_c(1, s) = c1
            Grid % faces_c(2, s) = -n_bnd_cells
            Grid % region % at_cell(-n_bnd_cells) = face_sect_bnd(n_face_sect)

          ! Neither c1 nor c2 are boundary cells
          else
            Grid % faces_c(1, s) = min(c1, c2)
            Grid % faces_c(2, s) = max(c1, c2)

          end if

        end do
      end if

    end if
  end do
  if(the_end) then
    print *, '# ERROR: Could not find face-based data in: ', file_name
    print *, '# This error is critical!  Exiting!'
    stop
  end if

  print '(a34,i9)', ' # Boundary cells from face data: ', n_bnd_cells

  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Checking faces                                           '
  print '(a60)', ' #----------------------------------------------------------'

  !----------------------------!
  !   Check for faces' nodes   !
  !----------------------------!
  do s = 1, Grid % n_faces
    do i_nod = 1, Grid % faces_n_nodes(s)
      n = Grid % faces_n(i_nod, s)
      if(n < 1) then
        print *, '# TROUBLE: some faces'' nodes have indexes less than 1!'
        print *, '# This error is critical.  Exiting now.!'
        stop
      end if
    end do
  end do

  !----------------------------------------!
  !   Check for duplicate nodes in faces   !
  !----------------------------------------!
  do s = 1, Grid % n_faces
    n = Grid % faces_n_nodes(s)
    do i_nod = 1, n
      do j_nod = i_nod+1, n
        if(Grid % faces_n(i_nod, s) .eq. Grid % faces_n(j_nod, s)) then
          print *, '# ERROR!  Duplicate nodes in face: ', s
          print *, '# This error is critical, exiting! '
          stop
        end if
      end do
    end do
  end do
  print *, '# No duplicate nodes in face data found, good!'

  !----------------------------------------!
  !   Check for duplicate cells in faces   !
  !----------------------------------------!
  do s = 1, Grid % n_faces
    if(Grid % faces_c(1, s) .eq. Grid % faces_c(2, s)) then
      print *, '# ERROR!  Duplicate cells in face: ', s
      print *, '# This error is critical, exiting! '
      stop
    end if
  end do
  print *, '# No duplicate cells in face data found, good!'

  !--------------------------------!
  !                                !
  !                                !
  !    Reconstruct cells' nodes    !
  !   (For non-polyhedral cells)   !
  !                                !
  !                                !
  !--------------------------------!
  print '(a60)', ' #=========================================================='
  print '(a60)', ' # Reconstructing cells (determining their nodes)           '
  print '(a60)', ' #----------------------------------------------------------'
  allocate(cell_visited_from(Grid % n_cells));
  cell_visited_from(:) = 0

  !---------------------------------------------!
  !   Handle all boundary cells to start with   !
  !---------------------------------------------!
  do s = 1, Grid % n_faces
    c2 = Grid % faces_c(2, s)

    ! It is a boundary cell, just copy nodes
    if(c2 .lt. 0) then
      n = Grid % faces_n_nodes(s)
      Grid % cells_n_nodes(c2) = n
      call Adjust_First_Dim(n, Grid % cells_n)
      Grid % cells_n(1:n, c2)  = Grid % faces_n(1:n, s)
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
  do s = 1, Grid % n_faces

    !--------------------------------------!
    !   Handle cells c1 and c2 in one go   !
    !--------------------------------------!
    do i_cel = 1, 2

      c = Grid % faces_c(i_cel, s)

      if(c .gt. 0) then

        !------------------------------------------------------------!
        !   First visit to pyramid or hexa from quadrilateral face   !
        !    For hexahedra that's face 5, for pyramid it's face 1    !
        !------------------------------------------------------------!
        if(Grid % faces_n_nodes(s) .eq. 4) then
          if( (Grid % cells_n_nodes(c) .eq. 5  .or.   &
               Grid % cells_n_nodes(c) .eq. 8) .and.  &
              cell_visited_from(c) .eq. 0) then
            Grid % cells_n(1:4, c) = Grid % faces_n(1:4, s)
            cell_visited_from(c) = s
          end if
        end if

        !--------------------------------------------------------------!
        !   First visit to tetrahedron or wedge from triangular face   !
        !             For both shapes that would be face 1             !
        !--------------------------------------------------------------!
        if(Grid % faces_n_nodes(s) .eq. 3) then
          if( (Grid % cells_n_nodes(c) .eq. 4  .or.   &
               Grid % cells_n_nodes(c) .eq. 6) .and.  &
              cell_visited_from(c) .eq. 0) then
            Grid % cells_n(1:3, c) = Grid % faces_n(1:3, s)
            cell_visited_from(c) = s
          end if
        end if  ! face is triangular

      end if  ! c > 0

    end do  ! through c1 and c2
  end do  ! through faces

  !------------------------------------------!
  !   Check if all cells have been visited   !
  !------------------------------------------!
  do c = 1, Grid % n_cells
    if(Grid % cells_n_nodes(c) .ge. 4 .and.  &
       Grid % cells_n_nodes(c) .le. 8) then
      if(cell_visited_from(c) .eq. 0) then
        print *, '# ERROR: Some cells have not been visited!'
        print *, '# This error is critical.  Exiting now.!'
        stop
      end if
    end if
  end do

  !---------------------------------------------------!
  !                                                   !
  !   Browse through all faces for the second visit   !
  !                                                   !
  !---------------------------------------------------!
  do s = 1, Grid % n_faces

    !--------------------------------------!
    !   Handle cells c1 and c2 in one go   !
    !--------------------------------------!
    do i_cel = 1, 2

      c = Grid % faces_c(i_cel, s)

      if(c .gt. 0) then

        !---------------------------!
        !   Face is quadrilateral   !
        !---------------------------!
        if(Grid % faces_n_nodes(s) .eq. 4) then

          !--------------------------------!
          !   Second visit to hexahedron   !
          !--------------------------------!
          if( Grid % cells_n_nodes(c) .eq. 8 .and.  &
              cell_visited_from(c) .ne.  s) then

            ! i_nod and j_nod, two consecutive nodes in the face
            do i_nod = 1, 4  ! nodes in face
              j_nod = i_nod + 1;  if(j_nod > 4) j_nod = j_nod - 4
              k_nod = i_nod + 2;  if(k_nod > 4) k_nod = k_nod - 4
              l_nod = i_nod + 3;  if(l_nod > 4) l_nod = l_nod - 4

              ! Face 1 same sense of rotation (see Cell_Numbering_Neu.f90)
              if(Grid % faces_n(i_nod, s) .eq. Grid % cells_n(1, c) .and.  &
                 Grid % faces_n(j_nod, s) .eq. Grid % cells_n(2, c) ) then
                Grid % cells_n(6, c) = Grid % faces_n(k_nod, s)
                Grid % cells_n(5, c) = Grid % faces_n(l_nod, s)
              end if

              ! Face 1 oposite sense of rotation (see Cell_Numbering_Neu.f90)
              if(Grid % faces_n(i_nod, s) .eq. Grid % cells_n(2, c) .and.  &
                 Grid % faces_n(j_nod, s) .eq. Grid % cells_n(1, c) ) then
                Grid % cells_n(5, c) = Grid % faces_n(k_nod, s)
                Grid % cells_n(6, c) = Grid % faces_n(l_nod, s)
              end if

              ! Face 3 same sense of rotation (see Cell_Numbering_Neu.f90)
              if(Grid % faces_n(i_nod, s) .eq. Grid % cells_n(3, c) .and.  &
                 Grid % faces_n(j_nod, s) .eq. Grid % cells_n(4, c) ) then
                Grid % cells_n(8, c) = Grid % faces_n(k_nod, s)
                Grid % cells_n(7, c) = Grid % faces_n(l_nod, s)
              end if

              ! Face 3 oposite sense of rotation (see Cell_Numbering_Neu.f90)
              if(Grid % faces_n(i_nod, s) .eq. Grid % cells_n(4, c) .and.  &
                 Grid % faces_n(j_nod, s) .eq. Grid % cells_n(3, c) ) then
                Grid % cells_n(7, c) = Grid % faces_n(k_nod, s)
                Grid % cells_n(8, c) = Grid % faces_n(l_nod, s)
              end if

            end do
          end if    ! if hexahedron

          !-----------------------------------------------------------!
          !   Second time you visit a wedge from quadrilateral face   !
          !-----------------------------------------------------------!
          if( Grid % cells_n_nodes(c) .eq. 6 .and.  &
              cell_visited_from(c) .ne.  s) then

            ! i_nod and j_nod, two consecutive nodes in the face
            do i_nod = 1, 4  ! nodes in face
              j_nod = i_nod + 1;  if(j_nod > 4) j_nod = j_nod - 4
              k_nod = i_nod + 2;  if(k_nod > 4) k_nod = k_nod - 4
              l_nod = i_nod + 3;  if(l_nod > 4) l_nod = l_nod - 4

              ! Face 1 same sense of rotation (see Cell_Numbering_Neu.f90)
              if(Grid % faces_n(i_nod, s) .eq. Grid % cells_n(1, c) .and.  &
                 Grid % faces_n(j_nod, s) .eq. Grid % cells_n(2, c) ) then
                Grid % cells_n(5, c) = Grid % faces_n(k_nod, s)
                Grid % cells_n(4, c) = Grid % faces_n(l_nod, s)
              end if

              ! Face 1 oposite sense of rotation (see Cell_Numbering_Neu.f90)
              if(Grid % faces_n(i_nod, s) .eq. Grid % cells_n(2, c) .and.  &
                 Grid % faces_n(j_nod, s) .eq. Grid % cells_n(1, c) ) then
                Grid % cells_n(4, c) = Grid % faces_n(k_nod, s)
                Grid % cells_n(5, c) = Grid % faces_n(l_nod, s)
              end if

              ! Face 2 same sense of rotation (see Cell_Numbering_Neu.f90)
              if(Grid % faces_n(i_nod, s) .eq. Grid % cells_n(2, c) .and.  &
                 Grid % faces_n(j_nod, s) .eq. Grid % cells_n(3, c) ) then
                Grid % cells_n(6, c) = Grid % faces_n(k_nod, s)
                Grid % cells_n(5, c) = Grid % faces_n(l_nod, s)
              end if

              ! Face 2 oposite sense of rotation (see Cell_Numbering_Neu.f90)
              if(Grid % faces_n(i_nod, s) .eq. Grid % cells_n(3, c) .and.  &
                 Grid % faces_n(j_nod, s) .eq. Grid % cells_n(2, c) ) then
                Grid % cells_n(5, c) = Grid % faces_n(k_nod, s)
                Grid % cells_n(6, c) = Grid % faces_n(l_nod, s)
              end if

            end do

          end if

        end if  ! face is quadrilateral

        !------------------------!
        !   Face is triangular   !
        !------------------------!
        if(Grid % faces_n_nodes(s) .eq. 3) then

          !-------------------------------------------------------------!
          !   Second time you visit a pyramid from triangular face      !
          !   (For pyramid, one node still missing, any face will do)   !
          !-------------------------------------------------------------!
          if( Grid % cells_n_nodes(c) .eq. 5 .and.  &
              cell_visited_from(c) .ne.  s   .and.  &
              cell_visited_from(c) .ne. -1) then
            do i_nod = 1, 3
              if(all(Grid % cells_n(1:4, c) .ne.  &
                     Grid % faces_n(i_nod, s))) then
                Grid % cells_n(5, c) = Grid % faces_n(i_nod, s)
              end if
            end do

            ! Don't visit this cell again
            cell_visited_from(c) = -1
          end if

          !-----------------------------------------------------------------!
          !   Second time you visit a tetrahedron from triangular face      !
          !   (For tetrahedron, one node still missing, any face will do)   !
          !-----------------------------------------------------------------!
          if( Grid % cells_n_nodes(c) .eq. 4 .and.  &
              cell_visited_from(c) .ne.  s   .and.  &
              cell_visited_from(c) .ne. -1) then
            do i_nod = 1, 3
              if(all(Grid % cells_n(1:3, c) .ne.  &
                     Grid % faces_n(i_nod, s))) then
                Grid % cells_n(4, c) = Grid % faces_n(i_nod, s)
              end if
            end do

            ! Don't visit this cell again
            cell_visited_from(c) = -1
          end if
        end if  ! face is triangular

      end if  ! c > 0

    end do  ! through c1 and c2
  end do  ! through s, faces

  !-----------------------!
  !                       !
  !   Find cells' faces   !
  !                       !
  !-----------------------!
  call Grid % Find_Cells_Faces()

  !--------------------------------------------------------------------!
  !                                                                    !
  !   With faces counted and stored, store nodes for each polyhedron   !
  !                                                                    !
  !--------------------------------------------------------------------!
  do c = 1, Grid % n_cells

    ! Only do this for polyhedral cells
    ! (For the other it was done above)
    if(Grid % cells_n_nodes(c) .eq. -1) then

      ! Accumulate nodes from all faces surrounding the cell
      n = 0
      do i_fac = 1, Grid % cells_n_faces(c)
        s = Grid % cells_f(i_fac, c)           ! take true face index
        do i_nod = 1, Grid % faces_n_nodes(s)
          n = n + 1
          all_nodes(n) = Grid % faces_n(i_nod, s)
        end do
      end do

      ! Perform a unique sort to remove duplicates
      call Sort % Unique_Int(all_nodes(1:n), n)

      call Adjust_First_Dim(n, Grid % cells_n)
      Grid % cells_n(1:n, c)  = all_nodes(1:n)
      Grid % cells_n_nodes(c) = -n

    end if  ! if cell was polyhedral

  end do

  !----------------------------!
  !   Check for cells' nodes   !
  !----------------------------!
  do c = 1, Grid % n_cells
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      if(n < 1) then
        print *, '# ERROR: some cells'' nodes have indexes less than 1!'
        print *, '# This error is critical.  Exiting now.!'
        stop
      end if
    end do
  end do

  !----------------------------------------!
  !   Check for duplicate nodes in cells   !
  !----------------------------------------!
  do c = 1, Grid % n_cells
    n = abs(Grid % cells_n_nodes(c))
    do i_nod = 1, n
      do j_nod = i_nod+1, n
        if(Grid % cells_n(i_nod, c) .eq. Grid % cells_n(j_nod, c)) then
          print *, '# ERROR!  Duplicate nodes in cell: ', c
          print *, '# This error is critical, exiting! '
          stop
        end if
      end do
    end do
  end do
  print *, '# No duplicate nodes in cell data found, good!'

  deallocate(cell_visited_from)

  !------------------------------------------------------------!
  !                                                            !
  !                                                            !
  !   Browse through all sections (fluid, boundary, and what   !
  !    ever) to store number of boundary conditions to Grid    !
  !                                                            !
  !                                                            !
  !------------------------------------------------------------!
  rewind(fu)
  the_end = .false.
  do while(.not. the_end)
    call File % Read_Line(fu, reached_end=the_end)
    if(.not. the_end .and. Line % n_tokens > 1) then
      if(Line % tokens(1) .eq. '(39' .or.  &
         Line % tokens(1) .eq. '(45') then

        ! Extract this face section position in the Fluent's mesh
        one_token = Line % tokens(2)
        l = len_trim(one_token)
        read(one_token(2:l), '(i16)') pos

        ! Extract boundary condition name
        if(index(Line % tokens(4), ')') .ne. 0) then  ! ')' is in the string
          one_token = Line % tokens(4)(1:index(Line % tokens(4), ')')-1)
        else
          one_token = Line % tokens(4)
        end if

        do n = 1, 2048
          if(face_sect_bnd(n) .ne. 0) then
            if(face_sect_pos(n) .eq. pos) then
              call String % To_Upper_Case(one_token)
              Grid % region % name(face_sect_bnd(n)) = one_token
            end if
          end if
        end do

      end if
    end if
  end do

  close(fu)

  call Profiler % Stop('Load_Fluent')

  end subroutine
