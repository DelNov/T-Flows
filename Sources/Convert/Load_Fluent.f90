!==============================================================================!
  subroutine Load_Fluent(grid)
!------------------------------------------------------------------------------!
!   Reads the Fluent's file format.                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  character(DL)        :: name_in
  integer              :: n_tri, n_quad, n_tetra, n_hexa, n_pyra, n_wedge
  integer              :: n_cells, n_bnd_cells, n_faces, n_nodes, n_face_nodes
  integer              :: c, c1, c2, c12, s, n, cell_type, fu, i, i_nod, j_nod
  integer              :: cell_s, cell_e, side_s, side_e, node_s, node_e
  integer, allocatable :: cell_visited_from(:)
  real                 :: dist(4)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MSH_MIXED = 0
  integer, parameter :: MSH_TRI   = 1
  integer, parameter :: MSH_TETRA = 2
  integer, parameter :: MSH_QUAD  = 3
  integer, parameter :: MSH_HEXA  = 4
  integer, parameter :: MSH_PYRA  = 5
  integer, parameter :: MSH_WEDGE = 6
!==============================================================================!

  call File_Mod_Set_Name(name_in, extension='.msh')
  call File_Mod_Open_File_For_Reading(name_in, fu)

  !--------------------------------------------------------!
  !                                                        !
  !                                                        !
  !   Read number of nodes, cells, faces, boundary cells   !
  !                                                        !
  !                                                        !
  !--------------------------------------------------------!
  print *, '# Reading header'

  !------------------------------------------!
  !   Look for number of nodes and read it   !
  !     Nodes' section starts with '(10'     !
  !------------------------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then
      if(line % tokens(1) .eq. '(10' .and. line % tokens(2) .eq. '(0') then
        read(line % tokens(4), '(z160)') grid % n_nodes
        print *, '# Number of nodes in header: ', grid % n_nodes
        goto 1
      end if
    end if
  end do
1 continue

  !------------------------------------------!
  !   Look for number of cells and read it   !
  !     Cells' section starts with '(12'     !
  !------------------------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then
      if(line % tokens(1) .eq. '(12' .and. line % tokens(2) .eq. '(0') then
        read(line % tokens(4), '(z160)') grid % n_cells
        print *, '# Number of cells in header: ', grid % n_cells
        goto 2
      end if
    end if
  end do
2 continue

  !------------------------------------------!
  !   Look for number of faces and read it   !
  !     Faces' section starts with '(13'     !
  !------------------------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then
      if(line % tokens(1) .eq. '(13' .and. line % tokens(2) .eq. '(0') then
        read(line % tokens(4), '(z160)') grid % n_faces
        print *, '# Number of faces in header: ', grid % n_faces
        goto 3
      end if
    end if
  end do
3 continue

  !------------------------------------------------!
  !   Read the faces to count the boundary cells   !
  !------------------------------------------------!
  print *, '# Reading face data to find number of boundary cells'
  n_bnd_cells = 0
  n_faces     = 0
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then

      ! Does the line mark the beginning of face-based data
      if(line % tokens(1) .eq. '(13' .and. line % tokens(2) .ne. '(0') then
        read(line % tokens(3), '(z160)') side_s  ! starting face
        read(line % tokens(4), '(z160)') side_e  ! ending face
        print *, '# Found a face zone from ', side_s, ' to ', side_e
        do s = side_s, side_e
          n_faces = n_faces + 1
          call File_Mod_Read_Line(fu)

          ! Number of nodes in this face
          read(line % tokens(1), *) n_face_nodes

          ! Read cells surrounding the face
          if     (n_face_nodes .eq. 3) then
            read(line % tokens(5), '(z160)') c1
            read(line % tokens(6), '(z160)') c2
          else if(n_face_nodes .eq. 4) then
            read(line % tokens(6), '(z160)') c1
            read(line % tokens(7), '(z160)') c2
          end if

          ! This was a boundary face
          if( (c1 .eq. 0) .or. (c2 .eq. 0) ) then
            n_bnd_cells = n_bnd_cells + 1
          end if

        end do
      end if

      ! Check if you read all the nodes
      if(n_faces .ge. grid % n_faces) goto 4

    end if
  end do
4 continue
  grid % n_bnd_cells = n_bnd_cells
  print *, '# Number of boundary cells from face data: ', grid % n_bnd_cells

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
  print *, '# Reading node coordinates'
  n_nodes = 0
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then

      ! Does the line mark the beginning of node-based data
      if(line % tokens(1) .eq. '(10' .and. line % tokens(2) .ne. '(0') then
        read(line % tokens(3), '(z160)') node_s  ! starting node
        read(line % tokens(4), '(z160)') node_e  ! ending node
        print *, '# Found a node zone from ', node_s, ' to ', node_e
        do n = node_s, node_e
          n_nodes = n_nodes + 1
          call File_Mod_Read_Line(fu)
          read(line % tokens(1), *)  grid % xn(n)
          read(line % tokens(2), *)  grid % yn(n)
          read(line % tokens(3), *)  grid % zn(n)
        end do
      end if

      ! Check if you read all the nodes
      if(n_nodes .ge. grid % n_nodes) goto 5

    end if
  end do
5 continue

  !-------------------------!
  !                         !
  !                         !
  !   Read the cell types   !  (Works only for one zone of mixed cell types!!!)
  !                         !
  !                         !
  !-------------------------!
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then
      if(line % tokens(1) .eq. '(12' .and. line % tokens(2) .ne. '(0') then

        read(line % tokens(3), '(z160)') cell_s  ! starting face
        read(line % tokens(4), '(z160)') cell_e  ! ending face
        print *, '# Found a cell zone from ', cell_s, ' to ', cell_e

        ! Initialize all counters
        n_cells = 0
        n_tri   = 0;  n_quad  = 0;  n_tetra = 0
        n_hexa  = 0;  n_pyra  = 0;  n_wedge = 0

6       continue

        call File_Mod_Read_Line(fu)

        do i = 1, line % n_tokens
          read(line % tokens(i), *) cell_type

          ! Update the counter for all cells
          n_cells = n_cells + 1

          ! Update the counters for cell types
          if(cell_type .eq. MSH_TRI) then
            n_tri = n_tri + 1
            grid % cells_n_nodes(n_cells) = 3
            grid % cells_n_faces(n_cells) = 1

          else if(cell_type .eq. MSH_QUAD) then
            n_quad = n_quad + 1
            grid % cells_n_nodes(n_cells) = 4
            grid % cells_n_faces(n_cells) = 1

          else if(cell_type .eq. MSH_TETRA) then
            n_tetra = n_tetra + 1
            grid % cells_n_nodes(n_cells) = 4
            grid % cells_n_faces(n_cells) = 4

          else if(cell_type .eq. MSH_HEXA)  then
            n_hexa = n_hexa + 1
            grid % cells_n_nodes(n_cells) = 8
            grid % cells_n_faces(n_cells) = 6

          else if(cell_type .eq. MSH_PYRA)  then
            n_pyra = n_pyra + 1
            grid % cells_n_nodes(n_cells) = 5
            grid % cells_n_faces(n_cells) = 5

          else if(cell_type .eq. MSH_WEDGE) then
            n_wedge = n_wedge + 1
            grid % cells_n_nodes(n_cells) = 6
            grid % cells_n_faces(n_cells) = 5

          else
            print *, '# ERROR: Unsupported cell type.'
            print *, '# This error is critical.  Exiting!'
            stop

          end if
        end do

        if(n_cells < grid % n_cells) then
          goto 6
        else
          goto 7
        end if

      end if
    end if
  end do
7 continue
  print *, '# Number of triangles:      ', n_tri
  print *, '# Number of quadrilaterals: ', n_quad
  print *, '# Number of tetrahedra:     ', n_tetra
  print *, '# Number of hexahedra:      ', n_hexa
  print *, '# Number of pyramids:       ', n_pyra
  print *, '# Number of wedges:         ', n_wedge

  !----------------------------------------------------!
  !                                                    !
  !                                                    !
  !   Read the faces again to store nodes, c1 and c2   !
  !                                                    !
  !                                                    !
  !----------------------------------------------------!
  print *, '# Reading face data to store c1, c2 and boundary cells'
  n_bnd_cells = 0
  n_faces     = 0
  rewind(fu)
  do
    call File_Mod_Read_Line(fu)
    if(line % n_tokens > 1) then

      ! Does the line mark the beginning of face-based data
      if(line % tokens(1) .eq. '(13' .and. line % tokens(2) .ne. '(0') then
        read(line % tokens(3), '(z160)') side_s  ! starting face
        read(line % tokens(4), '(z160)') side_e  ! ending face
        print *, '# Found a face zone from ', side_s, ' to ', side_e
        do s = side_s, side_e

          ! Increase the face counter
          n_faces = n_faces + 1
          call File_Mod_Read_Line(fu)

          ! Read and store number of nodes in this face
          read(line % tokens(1), *) n_face_nodes
          grid % faces_n_nodes(s) = n_face_nodes

          ! Read cells surrounding the face
          if     (n_face_nodes .eq. 3) then
            read(line % tokens(2), '(z160)') grid % faces_n(1, s)
            read(line % tokens(3), '(z160)') grid % faces_n(2, s)
            read(line % tokens(4), '(z160)') grid % faces_n(3, s)
            read(line % tokens(5), '(z160)') c1
            read(line % tokens(6), '(z160)') c2
          else if(n_face_nodes .eq. 4) then
            read(line % tokens(2), '(z160)') grid % faces_n(1, s)
            read(line % tokens(3), '(z160)') grid % faces_n(2, s)
            read(line % tokens(4), '(z160)') grid % faces_n(3, s)
            read(line % tokens(5), '(z160)') grid % faces_n(4, s)
            read(line % tokens(6), '(z160)') c1
            read(line % tokens(7), '(z160)') c2
          end if

          ! Case when c1 is a boundary cell
          if(c1 .eq. 0) then
            n_bnd_cells = n_bnd_cells + 1
            grid % faces_c(1, n_faces) = c2
            grid % faces_c(2, n_faces) = -n_bnd_cells

          ! Case when c2 is a boundary cell
          else if(c2 .eq. 0) then
            n_bnd_cells = n_bnd_cells + 1
            grid % faces_c(1, n_faces) = c1
            grid % faces_c(2, n_faces) = -n_bnd_cells

          ! Neither c1 nor c2 are boundary cells
          else
            grid % faces_c(1, n_faces) = min(c1, c2)
            grid % faces_c(2, n_faces) = max(c1, c2)

          end if

        end do
      end if

      ! Check if you read all the nodes
      if(n_faces .ge. grid % n_faces) goto 8

    end if
  end do
8 continue

  print *, '# Number of boundary cells from face data: ', n_bnd_cells

  !------------------------------!
  !                              !
  !                              !
  !   Reconstruct cells' nodes   !
  !                              !
  !                              !
  !------------------------------!
  print *, '# Reconstructing cells (determining their nodes)'
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
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    !---------------------------!
    !   Face is quadrilateral   !
    !---------------------------!
    if(grid % faces_n_nodes(s) .eq. 4) then

      !------------------------------------------------------------------------!
      !   First visit to pyramid or hexahedron from quadrilateral face as c1   !
      !------------------------------------------------------------------------!
      if( (grid % cells_n_nodes(c1) .eq. 5  .or.   &
           grid % cells_n_nodes(c1) .eq. 8) .and.  &
          cell_visited_from(c1) .eq. 0) then

        ! Check Cell_Numbering_Neu, there is a twist at 3, 4
        ! (Still not sure if that is the way it should be, not checked yet)
        grid % cells_n(1, c1) = grid % faces_n(1, s)
        grid % cells_n(2, c1) = grid % faces_n(2, s)
        grid % cells_n(4, c1) = grid % faces_n(3, s)
        grid % cells_n(3, c1) = grid % faces_n(4, s)
        cell_visited_from(c1) = s
      end if

      if(c2 .gt. 0) then

        !------------------------------------------------------------------!
        !   First visit to pyramid or hexa from quadrilateral face as c2   !
        !------------------------------------------------------------------!
        if( (grid % cells_n_nodes(c2) .eq. 5  .or.   &
             grid % cells_n_nodes(c2) .eq. 8) .and.  &
            cell_visited_from(c2) .eq. 0) then

          ! Check Cell_Numbering_Neu, there is a twist at 3, 4
          ! (Still not sure if that is the way it should be, not checked yet)
          grid % cells_n(1, c2) = grid % faces_n(1, s)
          grid % cells_n(2, c2) = grid % faces_n(2, s)
          grid % cells_n(4, c2) = grid % faces_n(3, s)
          grid % cells_n(3, c2) = grid % faces_n(4, s)
          cell_visited_from(c2) = s
        end if
      end if

    end if

    !------------------------!
    !   Face is triangular   !
    !------------------------!
    if(grid % faces_n_nodes(s) .eq. 3) then

      !--------------------------------------------------------------------!
      !   First visit to tetrahedron or wedge from triangular face as c1   !
      !--------------------------------------------------------------------!
      if( (grid % cells_n_nodes(c1) .eq. 4  .or.   &
           grid % cells_n_nodes(c1) .eq. 6) .and.  &
          cell_visited_from(c1) .eq. 0) then
        grid % cells_n(1:3, c1) = grid % faces_n(1:3, s)
        cell_visited_from(c1) = s

        ! Correct if first face is a boundary face
        ! (This was important for wedges and makes no difference for tets
        !  Should something like that be done also for quad faces, not sure)
        if(c2 < 0) then
          call Swap_Int(grid % cells_n(2, c1), grid % cells_n(3, c1))
        end if

      end if

      if(c2 .gt. 0) then

        !--------------------------------------------------------------------!
        !   First visit to tetrahedron or wedge from triangular face as c2   !
        !--------------------------------------------------------------------!
        if( (grid % cells_n_nodes(c2) .eq. 4  .or.   &
             grid % cells_n_nodes(c2) .eq. 6) .and.  &
            cell_visited_from(c2) .eq. 0) then
          grid % cells_n(1:3, c2) = grid % faces_n(1:3, s)
          cell_visited_from(c2) = s
        end if
      end if

    end if

  end do  ! through faces

  !---------------------------------------------------!
  !                                                   !
  !   Browse through all faces for the second visit   !
  !                                                   !
  !---------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    !---------------------------!
    !   Face is quadrilateral   !
    !---------------------------!
    if(grid % faces_n_nodes(s) .eq. 4) then

      !--------------------------------------!
      !   Second visit to hexahedron in c1   !
      !--------------------------------------!
      if( grid % cells_n_nodes(c1) .eq. 8 .and.  &
          cell_visited_from(c1) .ne.  s   .and.  &
          cell_visited_from(c1) .ne. -1) then

        ! Check if some nodes are already matching to make sure that
        ! you are at the opposite side of face which has been visited
        ! (Hence: skip this face if some nodes are already matching)
        do i_nod = 1, 4
          do j_nod = 1, 4
            if(grid % cells_n(i_nod, c1) .eq. grid % faces_n(j_nod, s)) then
              goto 10
            end if
          end do
        end do

        ! This is a wee-bit cumbersome, but couldn't think of anything else now
        ! (It associates stored nodes with new ones by their closeness)
        do i_nod = 1, 4  ! nodes already in hexahedral cell
          do j_nod = 1, 4  ! nodes in face
            dist(j_nod) = Math_Mod_Distance(                     &
                          grid % xn(grid % cells_n(i_nod, c1)),  &
                          grid % yn(grid % cells_n(i_nod, c1)),  &
                          grid % zn(grid % cells_n(i_nod, c1)),  &
                          grid % xn(grid % faces_n(j_nod,  s)),  &
                          grid % yn(grid % faces_n(j_nod,  s)),  &
                          grid % zn(grid % faces_n(j_nod,  s)))
          end do
          do j_nod = 1, 4
            if(minval(dist(1:4)) .eq. dist(j_nod)) then
              grid % cells_n(i_nod + 4, c1) = grid % faces_n(j_nod, s)
            end if
          end do

        end do

        !---------------------------------!
        !   Don't visit this cell again   !
        !---------------------------------!
        cell_visited_from(c1) = -1

      end if
10    continue  ! found matching nodes

      if(c2 .gt. 0) then

        !--------------------------------------!
        !   Second visit to hexahedron in c1   !
        !--------------------------------------!
        if( grid % cells_n_nodes(c2) .eq. 8 .and.  &
            cell_visited_from(c2) .ne.  s   .and.  &
            cell_visited_from(c2) .ne. -1) then

          ! Check if some nodes are already matching to make sure that
          ! you are at the opposite side of face which has been visited
          ! (Hence: skip this face if some nodes are already matching)
          do i_nod = 1, 4
            do j_nod = 1, 4
              if(grid % cells_n(i_nod, c2) .eq. grid % faces_n(j_nod, s)) then
                goto 20
              end if
            end do
          end do

          ! This is a wee-bit cumbersome, but couldn't think of anything else now
          ! (It associates stored nodes with new ones by their closeness)
          do i_nod = 1, 4    ! nodes already in hexahedral cell
            do j_nod = 1, 4  ! nodes in face
              dist(j_nod) = Math_Mod_Distance(                     &
                            grid % xn(grid % cells_n(i_nod, c2)),  &
                            grid % yn(grid % cells_n(i_nod, c2)),  &
                            grid % zn(grid % cells_n(i_nod, c2)),  &
                            grid % xn(grid % faces_n(j_nod,  s)),  &
                            grid % yn(grid % faces_n(j_nod,  s)),  &
                            grid % zn(grid % faces_n(j_nod,  s)))
            end do
            do j_nod = 1, 4
              if(minval(dist(1:4)) .eq. dist(j_nod)) then
                grid % cells_n(i_nod + 4, c2) = grid % faces_n(j_nod, s)
              end if
            end do

          end do

          !---------------------------------!
          !   Don't visit this cell again   !
          !---------------------------------!
          cell_visited_from(c2) = -1

        end if
20      continue  ! found matching nodes
      end if

    end if

    !------------------------!
    !   Face is triangular   !
    !------------------------!
    if(grid % faces_n_nodes(s) .eq. 3) then

      !------------------------------------------------------!
      !   Second visit to wedge from triangular face in c1   !
      !------------------------------------------------------!
      if( grid % cells_n_nodes(c1) .eq. 6 .and.  &
          cell_visited_from(c1) .ne.  s   .and.  &
          cell_visited_from(c1) .ne. -1) then

        ! This is a wee-bit cumbersome, but couldn't think of anything else now
        ! (It associates stored nodes with new ones by their closeness)
        do i_nod = 1, 3    ! nodes already in wedge cell
          do j_nod = 1, 3  ! nodes in face
            dist(j_nod) = Math_Mod_Distance(                     &
                          grid % xn(grid % cells_n(i_nod, c1)),  &
                          grid % yn(grid % cells_n(i_nod, c1)),  &
                          grid % zn(grid % cells_n(i_nod, c1)),  &
                          grid % xn(grid % faces_n(j_nod,  s)),  &
                          grid % yn(grid % faces_n(j_nod,  s)),  &
                          grid % zn(grid % faces_n(j_nod,  s)))
          end do
          do j_nod = 1, 3
            if(minval(dist(1:3)) .eq. dist(j_nod)) then
              grid % cells_n(i_nod + 3, c1) = grid % faces_n(j_nod, s)
            end if
          end do
        end do

        !---------------------------------!
        !   Don't visit this cell again   !
        !---------------------------------!
        cell_visited_from(c1) = -1

      end if

      if(c2 .gt. 0) then

        !--------------------------------------------------------------!
        !   Second time you visit a wedge from triangular face as c1   !
        !--------------------------------------------------------------!
        if( grid % cells_n_nodes(c2) .eq. 6 .and.  &
            cell_visited_from(c2) .ne.  s   .and.  &
            cell_visited_from(c2) .ne. -1) then

          ! This is a wee-bit cumbersome, but couldn't think of anything else now
          ! (It associates stored nodes with new ones by their closeness)
          do i_nod = 1, 3    ! nodes already in wedge cell
            do j_nod = 1, 3  ! nodes in face
              dist(j_nod) = Math_Mod_Distance(                     &
                            grid % xn(grid % cells_n(i_nod, c2)),  &
                            grid % yn(grid % cells_n(i_nod, c2)),  &
                            grid % zn(grid % cells_n(i_nod, c2)),  &
                            grid % xn(grid % faces_n(j_nod,  s)),  &
                            grid % yn(grid % faces_n(j_nod,  s)),  &
                            grid % zn(grid % faces_n(j_nod,  s)))
            end do
            do j_nod = 1, 3
              if(minval(dist(1:3)) .eq. dist(j_nod)) then
                grid % cells_n(i_nod + 3, c2) = grid % faces_n(j_nod, s)
              end if
            end do
          end do

          !---------------------------------!
          !   Don't visit this cell again   !
          !---------------------------------!
          cell_visited_from(c2) = -1

        end if
      end if

      !----------------------------------------------------------------!
      !   Second time you visit a pyramid from triangular face as c1   !
      !   (For pyramid, one node still missing, any face will do)      !
      !----------------------------------------------------------------!
      if( grid % cells_n_nodes(c1) .eq. 5 .and.  &
          cell_visited_from(c1) .ne.  s   .and.  &
          cell_visited_from(c1) .ne. -1) then
        do i_nod = 1, 3
          if(all(grid % cells_n(1:4, c1) .ne.  &
                 grid % faces_n(i_nod, s))) then
            grid % cells_n(5, c1) = grid % faces_n(i_nod, s)
          end if
        end do

        ! Don't visit this cell again
        cell_visited_from(c1) = -1
      end if

      if(c2 .gt. 0) then

        !----------------------------------------------------------------!
        !   Second time you visit a pyramid from triangular face as c2   !
        !   (For pyramid, one node still missing, any face will do)      !
        !----------------------------------------------------------------!
        if( grid % cells_n_nodes(c2) .eq. 5 .and.  &
            cell_visited_from(c2) .ne.  s   .and.  &
            cell_visited_from(c2) .ne. -1) then
          do i_nod = 1, 3
            if(all(grid % cells_n(1:4, c2) .ne.  &
                   grid % faces_n(i_nod, s))) then
              grid % cells_n(5, c2) = grid % faces_n(i_nod, s)
            end if
          end do

          ! Don't visit this cell again
          cell_visited_from(c2) = -1
        end if
      end if

      !--------------------------------------------------------------------!
      !   Second time you visit a tetrahedron from triangular face as c1   !
      !   (For tetrahedron, one node still missing, any face will do)      !
      !--------------------------------------------------------------------!
      if( grid % cells_n_nodes(c1) .eq. 4 .and.  &
          cell_visited_from(c1) .ne.  s   .and.  &
          cell_visited_from(c1) .ne. -1) then
        do i_nod = 1, 3
          if(all(grid % cells_n(1:3, c1) .ne.  &
                 grid % faces_n(i_nod, s))) then
            grid % cells_n(4, c1) = grid % faces_n(i_nod, s)
          end if
        end do

        ! Don't visit this cell again
        cell_visited_from(c1) = -1
      end if

      if(c2 .gt. 0) then

        !--------------------------------------------------------------------!
        !   Second time you visit a tetrahedron from triangular face as c2   !
        !   (For tetrahedron, one node still missing, any face will do)      !
        !--------------------------------------------------------------------!
        if( grid % cells_n_nodes(c2) .eq. 4 .and.  &
            cell_visited_from(c2) .ne.  s   .and.  &
            cell_visited_from(c2) .ne. -1) then
          do i_nod = 1, 3
            if(all(grid % cells_n(1:3, c2) .ne.  &
                   grid % faces_n(i_nod, s))) then
              grid % cells_n(4, c2) = grid % faces_n(i_nod, s)
            end if
          end do

          ! Don't visit this cell again
          cell_visited_from(c2) = -1
        end if
      end if

    end if

  end do

  deallocate(cell_visited_from)

  ! do c = 1, grid % n_cells
  !   print '(120i6)', grid % cells_n_nodes(c),  &
  !                    grid % cells_n(1:grid % cells_n_nodes(c), c)
  ! end do

  close(fu)

  end subroutine
