!==============================================================================!
  subroutine Print_Grid_Statistics(Grid, bounding_box)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)  :: Grid
  logical, optional :: bounding_box
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, s, n_cells, n_bnd_cells, n_faces, n_shadows
  integer :: n_tri, n_quad, n_tet, n_pyr, n_wed, n_hex, n_polg, n_polh
  integer :: max_n_polg, max_n_polh
  real    :: xmin, ymin, zmin, xmax, ymax, zmax  ! for bounding box
!==============================================================================!

  !------------!
  !   Header   !
  !------------!
  if(this_proc < 2) then
    print '(a)',  ' #========================================================='
    print '(a)',  ' #'
    print '(3a)', ' # Grid ', trim(Grid % name), ' statistics'
    print '(a)',  ' #'
    print '(a)',  ' #---------------------------------------------------------'
  end if

  !------------------!
  !   Bounding box   !
  !------------------!
  call Grid % Bounding_Box(xmin, ymin, zmin, xmax, ymax, zmax)

  if(this_proc < 2) then
    print '(a)',  ' # Bounding box:'
    print '(2(a,es10.3))', ' #   X range: ', xmin, '  to: ', xmax
    print '(2(a,es10.3))', ' #   Y range: ', ymin, '  to: ', ymax
    print '(2(a,es10.3))', ' #   Z range: ', zmin, '  to: ', zmax
  end if

  !-------------------------------------------------------------!
  !   Number of nodes, faces and cells, cells shape breakdown   !
  !-------------------------------------------------------------!
  ! Actually, I don't know how to work out number of nodes in parallel :-(

  n_faces    = 0
  n_shadows  = 0
  n_tri      = 0
  n_quad     = 0
  n_polg     = 0
  max_n_polg = 0  ! maximum number of nodes in polygon
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    if(Grid % comm % cell_proc(c1) .eq. this_proc) then
      n_faces = n_faces + 1
      max_n_polg = max(max_n_polg, Grid % faces_n_nodes(s))
      if(Grid % faces_s(s) .gt. s) then
        n_shadows = n_shadows + 1
      end if
      if(Grid % faces_n_nodes(s) .eq. 3) then
        n_tri = n_tri + 1
      else if(Grid % faces_n_nodes(s) .eq. 4) then
        n_quad = n_quad + 1
      else
        n_polg = n_polg + 1
      end if
    end if
  end do

  call Comm_Mod_Global_Sum_Int(n_faces)
  call Comm_Mod_Global_Sum_Int(n_shadows)
  call Comm_Mod_Global_Sum_Int(n_tri)
  call Comm_Mod_Global_Sum_Int(n_quad)
  call Comm_Mod_Global_Sum_Int(n_polg)
  call Comm_Mod_Global_Max_Int(max_n_polg)

  if(this_proc < 2) then
    print '(a)',     ' #- - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    print '(a,i9)',  ' # Number of faces:          ', n_faces
    print '(a,i9)',  ' # Number of shadows:        ', n_shadows
    print '(a)',     ' #'
    print '(a)',     ' # Face shape breakdown:'
    if(n_tri  > 0)  print '(a,i9)',  ' #   triangles:      ', n_tri
    if(n_quad > 0)  print '(a,i9)',  ' #   quadrilaterals: ', n_quad
    if(n_polg > 0)  print '(a,i9)',  ' #   polygons:       ', n_polg
    if(n_polg > 0) then
      print '(a)', ' #'
      print '(a,i3,a,i0.0,a)', ' # Maximum nodes in polygon:    ',  &
            max_n_polg, '  (of possible ', MAX_FACES_N_NODES, ')'
    end if
  end if

  n_bnd_cells = 0
  n_tri       = 0
  n_quad      = 0
  n_polg      = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % comm % cell_proc(c) .eq. this_proc) then
      n_bnd_cells = n_bnd_cells + 1
      if(Grid % cells_n_nodes(c) .eq. 3) then
        n_tri = n_tri + 1
      else if(Grid % cells_n_nodes(c) .eq. 4) then
        n_quad = n_quad + 1
      else
        n_polg = n_polg + 1
      end if
    end if
  end do

  n_cells = 0
  n_tet   = 0
  n_pyr   = 0
  n_wed   = 0
  n_hex   = 0
  n_polh  = 0
  max_n_polh = 0  ! maximum number of nodes in polyhedron
  do c = 1, Grid % n_cells
    if(Grid % comm % cell_proc(c) .eq. this_proc) then
      n_cells = n_cells + 1
      if(Grid % cells_n_nodes(c) .eq. 4) then
        n_tet = n_tet + 1
      else if(Grid % cells_n_nodes(c) .eq. 5) then
        n_pyr = n_pyr + 1
      else if(Grid % cells_n_nodes(c) .eq. 6) then
        n_wed = n_wed + 1
      else if(Grid % cells_n_nodes(c) .eq. 8) then
        n_hex = n_hex + 1
      else if(Grid % cells_n_nodes(c) < 0) then
        n_polh = n_polh + 1
        max_n_polh = max(max_n_polh, abs(Grid % cells_n_nodes(c))) ! needs abs!
      end if
    end if
  end do

  call Comm_Mod_Global_Sum_Int(n_bnd_cells)
  call Comm_Mod_Global_Sum_Int(n_cells)
  call Comm_Mod_Global_Sum_Int(n_tri)
  call Comm_Mod_Global_Sum_Int(n_quad)
  call Comm_Mod_Global_Sum_Int(n_polg)
  call Comm_Mod_Global_Sum_Int(n_tet)
  call Comm_Mod_Global_Sum_Int(n_pyr)
  call Comm_Mod_Global_Sum_Int(n_wed)
  call Comm_Mod_Global_Sum_Int(n_hex)
  call Comm_Mod_Global_Sum_Int(n_polh)
  call Comm_Mod_Global_Max_Int(max_n_polh)

  if(this_proc < 2) then
    print '(a)',     ' #- - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    print '(a,i9)',  ' # Number of bnd. cells:  ', n_bnd_cells
    print '(a,i9)',  ' # Number of inside cells:', n_cells
    print '(a)',     ' #'
    print '(a)',     ' # Cell shape breakdown:'
    if(n_tri  > 0)  print '(a,i9)',  ' #   triangles:      ', n_tri
    if(n_quad > 0)  print '(a,i9)',  ' #   quadrilaterals: ', n_quad
    if(n_polg > 0)  print '(a,i9)',  ' #   polygons:       ', n_polg
    if(n_tet  > 0)  print '(a,i9)',  ' #   tetrahedra:     ', n_tet
    if(n_pyr  > 0)  print '(a,i9)',  ' #   pyramids:       ', n_pyr
    if(n_wed  > 0)  print '(a,i9)',  ' #   wedges:         ', n_wed
    if(n_hex  > 0)  print '(a,i9)',  ' #   hexahedra:      ', n_hex
    if(n_polh > 0)  print '(a,i9)',  ' #   polyhedra:      ', n_polh

    if(n_polh > 0) then
      print '(a)', ' #'
      print '(a,i3,a,i0.0,a)', ' # Maximum nodes in polyhedron: ',  &
            max_n_polh, '  (of possible ', MAX_CELLS_N_NODES, ')'
    end if
  end if

  if(this_proc < 2) then
    print *, '#---------------------------------------------------------'
  end if

  end subroutine
