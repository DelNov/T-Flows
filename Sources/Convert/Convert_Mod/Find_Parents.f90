!==============================================================================!
  subroutine Find_Parents(Convert, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is critical for meshes where the connection between
!>  boundary cells and their corresponding internal (parent) cells is not
!>  directly given.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: Allocates memory for various arrays and initializes      !
!     variables. it starts with identifying all nodes on the boundary.         !
!   * Cells near boundaries: Counts the cells near the boundaries and          !
!     allocates working arrays for further processing.                         !
!   * Real work begins: This part involves the core logic of the subroutine.   !
!     It iterates over boundary regions and cells near boundaries to identify  !
!     potential parent-child relationships between boundary and internal cells !
!   * Sorting and matching process: The nodes of each potential face are       !
!     sorted and compared across boundary cells and internal cells to find     !
!     matches, which indicate a parent-child relationship.                     !
!   * Assigning parent information: Once a match is found, the subroutine      !
!     assigns the corresponding internal cell as the parent of the boundary    !
!     cell, updating the Grid % cells_bnd_region data structure.               !
!   * Final checks and warnings: The subroutine checks if the total number of  !
!     matches equals the number of boundary cells. if not, it issues a         !
!     warning, as this might indicate inconsistencies in the mesh.             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
!-----------------------------------[Locals]-----------------------------------!
  integer              :: fn(6,4), i_fac, j_fac, i_nod, i_cel, c1, c2
  integer              :: nodes(4), n_match, n_match_tot, dir, bc, cnt_c, cnt_f
  integer              :: n_face_nodes ! number of nodes in a face
  integer              :: n_cell_faces ! number of faces in a cell
  logical, allocatable :: is_node_bnd(:)
  integer, allocatable :: cell_near_bnd(:)
  integer, allocatable :: cr1(:), cr2(:), cr3(:), w1(:), w2(:)  ! for sorting
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  call Profiler % Start('Find_Parents')

  print *, '#================================================='
  print *, '# Parent information not given in the input file!'
  print *, '# Looking for parents. This may take a few minutes'
  print *, '#-------------------------------------------------'

  print *, '# Number of boundary cells: ', Grid % n_bnd_cells
  print *, '# Number of inside cells:   ', Grid % n_cells

  if(Grid % n_cells .eq. 0) then
    call Message % Error(50,                                &
      "Number of cells inside the domain is zero. \n "  //  &
      "Are you sure you meshed the domain in 3D?",          &
      file=__FILE__, line=__LINE__)
  end if

  allocate(is_node_bnd  (Grid % n_nodes));  is_node_bnd(:)   = .false.
  allocate(cell_near_bnd(Grid % n_cells));  cell_near_bnd(:) = 0

  !--------------------------------------------------------------!
  !   Find all cells near boundary and allocate working arrays   !
  !--------------------------------------------------------------!

  ! Mark all nodes on the boundary
  is_node_bnd(:) = .false.
  do c2 = -Grid % n_bnd_cells, -1
    if( Grid % region % at_cell(c2) .gt. 0 ) then
      do i_nod = 1, Grid % cells_n_nodes(c2)  ! 3 or 4
        is_node_bnd( Grid % cells_n(i_nod, c2) ) = .true.
      end do
    end if
  end do

  ! Count cells near boundaries
  cnt_c = 0
  cnt_f = Grid % n_bnd_cells
  do c1 = 1, Grid % n_cells
    do i_nod = 1, Grid % cells_n_nodes(c1)  ! 4 to 8
      if( is_node_bnd(Grid % cells_n(i_nod, c1)) ) then
        cnt_c = cnt_c + 1
        cell_near_bnd(cnt_c) = c1
        if(Grid % cells_n_nodes(c1) .eq. 4) cnt_f = cnt_f + 4  ! TET
        if(Grid % cells_n_nodes(c1) .eq. 5) cnt_f = cnt_f + 5  ! PYR
        if(Grid % cells_n_nodes(c1) .eq. 6) cnt_f = cnt_f + 5  ! WED
        if(Grid % cells_n_nodes(c1) .eq. 8) cnt_f = cnt_f + 6  ! HEX
        goto 1
      end if
    end do
1   continue
  end do
  print *, '# Number of all cells near boundary:          ', cnt_c
  print *, '# Number of all cells''s faces near boundary: ', cnt_f

  allocate(cr1(cnt_f));  cr1(:) = 0
  allocate(cr2(cnt_f));  cr2(:) = 0
  allocate(cr3(cnt_f));  cr3(:) = 0
  allocate(w1 (cnt_f));  w1 (:) = 0
  allocate(w2 (cnt_f));  w2 (:) = 0

  !----------------------!
  !   Real work begins   !
  !----------------------!
  n_match_tot = 0

  do bc = Boundary_Regions()

    cr1(:) = HUGE_INT
    cr2(:) = HUGE_INT
    cr3(:) = HUGE_INT
    w1 (:) = HUGE_INT
    w2 (:) = HUGE_INT

    cnt_f = 0
    do c2 = -Grid % n_bnd_cells, -1
      if( Grid % region % at_cell(c2) .eq. bc ) then

        ! Increase total face count
        cnt_f = cnt_f + 1

        ! Form and store local nodes list
        nodes(:) = HUGE_INT
        do i_nod = 1, Grid % cells_n_nodes(c2)  ! 3 or 4
          nodes(i_nod) = Grid % cells_n(i_nod, c2)
        end do
        call Sort % Int_Array(nodes(1:4))
        cr1(cnt_f) = nodes(1)
        cr2(cnt_f) = nodes(2)
        cr3(cnt_f) = nodes(3)
        w1 (cnt_f) = c2        ! not sure this will be needed, but OK
        w2 (cnt_f) = bc        ! not sure this will be needed, but OK
      end if
    end do

    do i_cel = 1, cnt_c
      c1 = cell_near_bnd(i_cel)

      ! Fetch the faces from this cell
      if(Grid % cells_n_nodes(c1) .eq. 4) fn = TET
      if(Grid % cells_n_nodes(c1) .eq. 5) fn = PYR
      if(Grid % cells_n_nodes(c1) .eq. 6) fn = WED
      if(Grid % cells_n_nodes(c1) .eq. 8) fn = HEX

      n_cell_faces = 6  ! assume it is a hexahedron
      if(Grid % cells_n_nodes(c1) .eq. 4) n_cell_faces = 4
      if(Grid % cells_n_nodes(c1) .eq. 5) n_cell_faces = 5
      if(Grid % cells_n_nodes(c1) .eq. 6) n_cell_faces = 5

      ! Browse through all possible faces and their nodes
      do i_fac = 1, n_cell_faces

        ! Increase total face count
        cnt_f = cnt_f + 1

        ! Estimate number of nodes for this face
        n_face_nodes = 4             ! assume it is a quad
        if( fn(i_fac, 4) < 0 ) then  ! nope, it was a triangle
          n_face_nodes = 3
        end if

        ! Form and store local nodes list
        nodes(:) = HUGE_INT
        do i_nod = 1, n_face_nodes
          if(fn(i_fac, i_nod) > 0) then  ! if this node exists
            nodes(i_nod) = Grid % cells_n(fn(i_fac, i_nod), c1)
          end if
        end do
        call Sort % Int_Array(nodes(1:4))
        cr1(cnt_f) = nodes(1)
        cr2(cnt_f) = nodes(2)
        cr3(cnt_f) = nodes(3)
        w1 (cnt_f) = c1
        w2 (cnt_f) = i_fac

      end do  ! i_fac

    end do  ! c1

    call Sort % Three_Int_Carry_Two_Int(                   &
                cr1(1:cnt_f), cr2(1:cnt_f), cr3(1:cnt_f),  &
                w1(1:cnt_f), w2(1:cnt_f))

    n_match = 0
    do i_fac = 1, cnt_f-1
      j_fac = i_fac + 1
      if(  sign(1, w1(i_fac))  &
         + sign(1, w1(j_fac)) .eq. 0 ) then  ! they have different signs

        ! Is it a match?
        if( cr1(i_fac) .eq. cr1(j_fac) .and.  &
            cr2(i_fac) .eq. cr2(j_fac) .and.  &
            cr3(i_fac) .eq. cr3(j_fac) ) then

          n_match = n_match + 1
          if( w1(i_fac) > 0 ) then  ! i_fac stores internal cell
            c1  = w1(i_fac)
            dir = w2(i_fac)
          else                      ! j_fac stores internal cell
            c1  = w1(j_fac)
            dir = w2(j_fac)
          end if
          Grid % cells_bnd_region(dir, c1) = bc

        end if
      end if
    end do

    print *, '# Number of matches: ', n_match
    n_match_tot = n_match_tot + n_match

  end do  ! bc

  if(n_match_tot .ne. Grid % n_bnd_cells) then
    call Message % Warning(80,                                                 &
              'Number of boundary cells found here is not the same as   '  //  &
              'the number of boundary cells prescribed in Gmsh. Possible ' //  &
              'cause for it is that some internal faces in Gmsh are '      //  &
              'marked as boundary conditions (physical groups). \n \n '    //  &
              'It might also be that you deleted some volumes in Gmsh '    //  &
              'intentionally to ensure conformal mappings for simulations '//  &
              'in multiple domains.  If that is the case, it is probably ' //  &
              'safe to ignore this warning.  \n \n '                       //  &
              'It is, nonetheless, advised to check the Gmsh mesh for '    //  &
              'physical groups.',                                              &
              file=__FILE__, line=__LINE__)
  end if

  call Profiler % Stop('Find_Parents')

  end subroutine
