!==============================================================================!
  subroutine Save_Fluent(Convert, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to save grid files in Fluent's file format.    !
!>  It uses the data from the provided Grid object to create into a .cas file, !
!>  which is compatible with the Fluent software.  For the time being, it      !
!>  saves the .cas file in the ASCII mode only.                                !
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * File opening and initial setup: Opens a new file in ASCII mode and       !
!     prepares it for writing grid data structured for Fluent compatibility.   !
!   * Header information writing: Writes basic header information to indicate  !
!     the exporting software (T-Flows) and set up the overall structure of the !
!     grid file including declarations of dimensions, nodes, cells, and faces. !
!   * Writing node data: Outputs all node coordinates in the grid to the file, !
!     recording each node's position (x, y, z).                                !
!   * Writing cell data: Outputs the type and count of each cell in the grid   !
!     according to Fluent's predefined mixed zone type.                        !
!   * Writing face data: Outputs the face data divided by boundary and inside  !
!     regions.  Each region is saved as Fluent's mixed zone and to each        !
!     reagion it assignes its unique number.                                   !
!   * For each face it writes the number of its nodes (needed for mixed        !
!     regions), the node indices and the indices of cells around the face.     !
!   * Final Processing: Closes the file and ends the write process.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert    !! parent class
  type(Grid_Type)     :: Grid       !! grid being saved
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
  character(SL) :: file_name
  integer       :: fu, n, c, i_nod, s, c1, c2, reg, reg_type, reg_zone
!==============================================================================!

  call Profiler % Start('Save_Fluent')

  !------------------------------------------!
  !                                          !
  !   Set the name of the file and open it   !
  !                                          !
  !------------------------------------------!
  file_name = trim(Grid % name)
  call String % To_Lower_Case(file_name)
  call File % Set_Name(file_name, extension=".cas")
  call File % Open_For_Writing_Ascii(file_name, fu)

  !------------------!
  !                  !
  !   Write header   !
  !                  !
  !------------------!

  !--------------------------------------------!
  !   Info on the software creating the file   !
  !--------------------------------------------!
  write(fu, '(a)') "(0 +---------------------------+)"
  write(fu, '(a)') "(0 |                           |)"
  write(fu, '(a)') "(0 |   Exported from T-Flows   |)"
  write(fu, '(a)') "(0 |                           |)"
  write(fu, '(a)') "(0 +---------------------------+)"

  write(fu, '(a)') "(0 File structure for the face-based grid"
  write(fu, '(a)') ""
  write(fu, '(a)') "   Declaration section:"
  write(fu, '(a)') "     Dimensions:  (2 3)"
  write(fu, '(a)') "     Nodes:       (10 (0 1 n_nodes 0))"
  write(fu, '(a)') "     Cells:       (12 (0 1 n_cells 0))"
  write(fu, '(a)') "     Faces:       (13 (0 1 n_faces 0))"
  write(fu, '(a)') ""
  write(fu, '(a)') "   Data sections:"
  write(fu, '(a)') "     Nodes:       (10 (zone_id 1 n_nodes type)"
  write(fu, '(a)') "                   x y z"
  write(fu, '(a)') "                   x y z"
  write(fu, '(a)') "                   ..."
  write(fu, '(a)') "                  )"
  write(fu, '(a)') "     Cells:       (12 (zone_id 1 n_cells type zone_type)"
  write(fu, '(a)') "                   cell_type_1 cell_type_2 ..."
  write(fu, '(a)') "                  )"
  write(fu, '(a)') "     Faces:       (13 (zone_id start end type zone_type)"
  write(fu, '(a)') "                   n v1 ... vn c1 c2"
  write(fu, '(a)') "                   n v1 ... vn c1 c2"
  write(fu, '(a)') "                   ..."
  write(fu, '(a)') "                  )"
  write(fu, '(a)') "   The zones are assigned as follows:"
  write(fu, '(a)') "     - Nodes are all in zone 1"
  write(fu, '(a)') "     - Cells are also in zone 1 and are of type 0 (mixed)"
  write(fu, '(a)') "     - Faces inside the domain are in zone 2"
  write(fu, '(a)') "     - Faces at boundaries are in zones 2 + boundary region"
  write(fu, '(a)') "     - All face zones are of type 0 (mixed type like cells)"
  write(fu, '(a)') "   Variable 'type' in all data sections is equal to 1"
  write(fu, '(a)') ")"

  !------------------------------------------------!
  !   Dimensions; it is a three-dimensional grid   !
  !------------------------------------------------!
  write(fu, '(a)') ""
  write(fu, '(a)') "(0 +------------------------------------+)"
  write(fu, '(a)') "(0 |   It is a three-dimensional grid   |)"
  write(fu, '(a)') "(0 +------------------------------------+)"
  write(fu, '(a)') "(2 3)"

  !--------------------------------------------!
  !   Write number of nodes, cells and faces   !
  !--------------------------------------------!
  write(fu, '(a)') ""
  write(fu, '(a)') "(0 +--------------------------------------+)"
  write(fu, '(a)') "(0 |   Number of nodes, cells and faces   |)"
  write(fu, '(a)') "(0 +--------------------------------------+)"
  write(fu, '(a, z9, a)')  "(10 (0 1 ", Grid % n_nodes, " 0))"
  write(fu, '(a, z9, a)')  "(12 (0 1 ", Grid % n_cells, " 0))"
  write(fu, '(a, z9, a)')  "(13 (0 1 ", Grid % n_faces, " 0))"

  !-----------------!
  !                 !
  !   Write nodes   !
  !                 !
  !-----------------!
  write(fu, '(a)') ""
  write(fu, '(a)') "(0 +------------------------+)"
  write(fu, '(a)') "(0 |                        |)"
  write(fu, '(a)') "(0 |   Nodes' coordinates   |)"
  write(fu, '(a)') "(0 |                        |)"
  write(fu, '(a)') "(0 +------------------------+)"
  write(fu, '(a, z9, a)')  "(10 (1 1 ", Grid % n_nodes, " 1) ("
  do n = 1, Grid % n_nodes
    write(fu, '(3es15.6)')  Grid % xn(n), Grid % yn(n), Grid % zn(n)
  end do
  write(fu, '(a)')  "))"

  !-----------------!
  !                 !
  !   Write cells   !
  !                 !
  !-----------------!
  write(fu, '(a)') ""
  write(fu, '(a)') "(0 +------------------+)"
  write(fu, '(a)') "(0 |                  |)"
  write(fu, '(a)') "(0 |   Cells' types   |)"
  write(fu, '(a)') "(0 |                  |)"
  write(fu, '(a)') "(0 +------------------+)"
  write(fu, '(a, 3z9, a)') "(12 (1 1 ", Grid % n_cells, 1, MIXED_ZONE, ") ("

  do c = 1, Grid % n_cells

    ! Take the number of nodes in the cell ...
    n = Grid % cells_n_nodes(c)

    ! ... and act accordingly
    if(n .lt. 0) then
      write(fu, '(z2)', advance="no")  CELL_POLY
    else if(n .eq. 4) then
      write(fu, '(z2)', advance="no")  CELL_TETRA
    else if(n .eq. 5) then
      write(fu, '(z2)', advance="no")  CELL_PYRA
    else if(n .eq. 6) then
      write(fu, '(z2)', advance="no")  CELL_WEDGE
    else if(n .eq. 8) then
      write(fu, '(z2)', advance="no")  CELL_HEXA
    else
      call Message % Error(80,                                             &
                           'Unknown cell type and I don''t know why. '//   &
                           'This error is critical. Exiting!',             &
                           file=__FILE__, line=__LINE__, one_proc=.true.)
    end if
  end do
  write(fu, '(a)')  "))"

  !-----------------!
  !                 !
  !   Write faces   !
  !                 !
  !-----------------!
  write(fu, '(a)') ""
  write(fu, '(a)') "(0 +----------------------------------------------+)"
  write(fu, '(a)') "(0 |                                              |)"
  write(fu, '(a)') "(0 |   Faces' nodes and cells, region by region   |)"
  write(fu, '(a)') "(0 |                                              |)"
  write(fu, '(a)') "(0 +----------------------------------------------+)"

  do reg = 1, Grid % n_bnd_regions + 1

    ! Specify region type
    reg_type = 3                              ! default is wall
    if(reg .eq. Grid % n_bnd_regions + 1)  &  ! correct for interior
      reg_type = 2

    ! Specify region zone
    reg_zone = reg + 2                        ! default are boundary regions
    if(reg .eq. Grid % n_bnd_regions + 1)  &  ! correct for interiors
      reg_zone = 2

    write(fu, '(a)')        ""
    write(fu, '(a)')        "(0 +----------------+)"
    write(fu, '(a, i3, a)') "(0 |   Region ", reg_zone, "   |)"
    if(reg_type .eq. 2) then
    write(fu, '(a, i3, a)') "(0 |   (interior)   |)"
    else
    write(fu, '(a, i3, a)') "(0 |   (boundary)   |)"
    end if
    write(fu, '(a)')        "(0 +----------------+)"
    ! Leading line of the new region
    write(fu, '(a, z0, z9, z9, z2, z2, a)')   &
               "(13 (",                       &
               reg_zone,                      &
               Grid % region % f_face(reg),   &
               Grid % region % l_face(reg),   &
               reg_type,                      &
               MIXED_ZONE,                    &
               ") ("

    !------------------------------------------!
    !   Browse through faces of that regions   !
    !------------------------------------------!
    do s = Grid % region % f_face(reg), Grid % region % l_face(reg)

      ! Write the number of nodes for this face
      write(fu, '(z9, 1x)', advance="no")  Grid % faces_n_nodes(c)

      ! Write all nodes for the face
      do i_nod = 1, Grid % faces_n_nodes(s)
        n = Grid % faces_n(i_nod, s)
        write(fu, '(z9, 1x)', advance="no")  n
      end do

      ! Write cells surroundig the face
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)
      c2 = max(0, c2)            ! if it is a boundary cell, set to zero

      write(fu, '(2z9)', advance="no")  c1, c2

      ! End the line for face s
      write(fu, '()')
    end do
    write(fu, '(a)')  "))"
  end do

  !--------------------!
  !                    !
  !   Close the file   !
  !                    !
  !--------------------!
  close(fu)

  call Profiler % Stop('Save_Fluent')

  end subroutine
