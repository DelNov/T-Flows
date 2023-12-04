!==============================================================================!
  subroutine Load_Obj(Convert, Grid, file_name)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to read and process a 3D object in the
!>  Wavefront OBJ file format, populating the Grid structure with its
!>  geometric information.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Opening the OBJ file: The subroutine starts by opening the specified     !
!     OBJ file in ASCII mode for reading.                                      !
!   * Counting nodes and faces:                                                !
!     - It counts the number of vertices (v) and faces (f) present in the      !
!       OBJ file.                                                              !
!     - These counts are used to allocate memory for storing node and face     !
!       information in the Grid structure.                                     !
!   * Reading coordinates and face nodes:                                      !
!     - The subroutine reads through the file again, this time extracting the  !
!       vertex coordinates and the indices of vertices that make up each face. !
!     - It swaps the y and z coordinates during this process.                  !
!   * Normalizing the Object:                                                  !
!     - It calculates the minimum and maximum values for the x, y, and z       !
!       coordinates of the vertices.                                           !
!     - The object is scaled based on its height (z dimension), and its origin !
!       is shifted to the center of its bounding box in the x and y directions ! 
!       This normalization ensures that the object fits within a standard      !
!       coordinate system and size range.                                      !
!   * Storing Geometry in Grid:                                                !
!     - The Grid structure is populated with the normalized node coordinates   !
!       and the faces' vertices.                                               !
!     - This information is used to represent the 3D geometry of the object    !
!       within T-Flows.                                                        !
!   * Handling of File End: The subroutine includes checks for the end of the  !
!     file to ensure that all necessary data is read without errors.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert    !! parent class
  type(Grid_Type)     :: Grid       !! grid being converted
  character(SL)       :: file_name  !! file name
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: tok
  integer       :: fu, s, n, i, i_nod
  logical       :: the_end
  real          :: t_min_x, t_max_x, t_min_y, t_max_y, t_min_z, t_max_z
  real          :: x_o, y_o
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  call Profiler % Start('Load_Obj')

  !------------------------!
  !                        !
  !   Read a single tree   !
  !                        !
  !------------------------!
  call File % Open_For_Reading_Ascii(file_name, fu)

  !---------------------------!
  !   Count nodes and faces   !
  !---------------------------!
  Grid % n_nodes = 0
  Grid % n_faces = 0
  Grid % name    = trim("trees")
  do
    call File % Read_Line(fu, reached_end=the_end)
    if(Line % tokens(1) == 'v') then
      Grid % n_nodes = Grid % n_nodes + 1
    end if
    if(Line % tokens(1) == 'f') then
      Grid % n_faces = Grid % n_faces + 1
    end if
    ! if(Line % tokens(1) .eq. '$MeshFormat') exit
    if(the_end) then
      goto 1
    end if
  end do
1 continue
  print '(a,i9)', ' # Nodes found: ', Grid % n_nodes
  print '(a,i9)', ' # Faces found: ', Grid % n_faces

  call Grid % Allocate_Nodes(Grid % n_nodes)
  call Grid % Allocate_Faces(Grid % n_faces, 0)

  rewind(fu)

  !---------------------------------------!
  !   Read coordinates and faces' nodes   !
  !---------------------------------------!
  n = 0
  s = 0
  do
    call File % Read_Line(fu, reached_end=the_end)
    if(Line % tokens(1) == 'v') then
      n = n + 1
      ! Swap y and z
      read(Line % tokens(2), *) Grid % xn(n)
      read(Line % tokens(3), *) Grid % zn(n)
      read(Line % tokens(4), *) Grid % yn(n)
    end if
    if(Line % tokens(1) == 'f') then
      s = s + 1
      Grid % faces_n_nodes(s) = Line % n_tokens - 1
      call Enlarge % Matrix_Int(Grid % faces_n, i=Grid % faces_n_nodes(s))
      do i_nod = 1, Grid % faces_n_nodes(s)
        tok = Line % tokens(i_nod + 1)
        i = index(tok, '/')
        read(tok(1:i), *) Grid % faces_n(i_nod, s)
      end do
    end if
    if(the_end) then
      goto 2
    end if
  end do
2 continue

  close(fu)

  !-------------------------------!
  !   Normalize the single tree   !
  !-------------------------------!
  t_min_x = minval(Grid % xn)
  t_max_x = maxval(Grid % xn)
  t_min_y = minval(Grid % yn)
  t_max_y = maxval(Grid % yn)
  t_min_z = minval(Grid % zn)
  t_max_z = maxval(Grid % zn)

  ! Scale the tree by height
  Grid % zn(:) = Grid % zn - t_min_z
  Grid % xn(:) = Grid % xn / (t_max_z - t_min_z)
  Grid % yn(:) = Grid % yn / (t_max_z - t_min_z)
  Grid % zn(:) = Grid % zn / (t_max_z - t_min_z)

  t_min_x = minval(Grid % xn)
  t_max_x = maxval(Grid % xn)
  t_min_y = minval(Grid % yn)
  t_max_y = maxval(Grid % yn)
  x_o = 0.5 * (t_max_x + t_min_x)
  y_o = 0.5 * (t_max_y + t_min_y)
  Grid % xn(:) = Grid % xn(:) - x_o
  Grid % yn(:) = Grid % yn(:) - y_o

  call Profiler % Stop('Load_Obj')

  end subroutine
