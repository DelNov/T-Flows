!==============================================================================!
  subroutine Insert_Buildings(Convert, Grid)
!------------------------------------------------------------------------------!
!>  This soubroutine was developed to import mesh with buildings, and place
!>  them both on a ground defined with an STL file.  This subroutine was
!>  traditionally not a part of the standard T-Flows distribution, but after
!>  publishing a couple of papers with it, we decided to include it.  It is
!>  rather long and often chaning, and it will probably stay like that in the
!>  foreseeable future.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_N          = 64
  integer, parameter :: ITERS          = 24
  integer, parameter :: FOR_BUILDINGS  =  1
  integer, parameter :: FOR_CHIMNEYS   =  2
  integer, parameter :: FOR_POROSITIES =  3
!-----------------------------------[Locals]-----------------------------------!

  ! For initial check
  logical :: has_top

  ! Variables for phase I: reading the ground definition from an STL file
  character(len=SL) :: name_in
  integer           :: n_ground_facets
  real, allocatable :: ground_x(:,:), ground_y(:,:), ground_z(:,:)
  integer           :: f, fu, i_ver
  real              :: x_min, x_max, y_min, y_max, z_min, z_max, area

  ! New variables for phase II: sorting nodes and cells, align where possible
  integer :: a_node_layers, n_node_layers  ! absolute and naive layers
  integer :: a_cell_layers, n_cell_layers  ! absolute and naive layers
  integer :: c, n, i_nod, nu

  ! New variables for phase III: searching for cells in buildings
  ! (It's more accurate doing this before placing grid on ground)
  integer, allocatable :: cell_in_building(:), node_on_building(:)
  integer, allocatable :: cell_in_chimney(:),  node_on_chimney(:)
  character(len=SL)    :: bc_name
  real                 :: height                  ! building height
  integer              :: run, region, dir, cu, cnt, bc
  integer              :: n_buildings
  integer              :: n_porosities            ! for porosities

  ! New variables for phase IV: place grid on the ground defined in STL
  real :: area_1, area_2, area_3, z_ground, z_scale, z_from_top

  ! New variables for phase V: flatten the tops of the buildings
  integer, allocatable :: buildings_n_ground_nodes(:)
  integer, allocatable :: buildings_ground_nodes(:,:)
  integer              :: b, lay, top, iter
  real                 :: z_avg

  ! New variables for phase VI: manage boundary condition names and numbers
  integer :: ground_bc

  ! New variables for phase VII: store only cells which are not in buildings
  integer, allocatable :: i_work_1(:)
  integer, allocatable :: i_work_2(:,:)
  integer, allocatable :: i_work_3(:,:)
  integer, allocatable :: i_work_4(:)

  ! New variables for phase VIII: insert new boundary faces
  integer :: fn(6,4), n_f_nod, f_nod(4)

  ! New variables for phase IX: smooth the grid lines in z direction
  integer, allocatable :: nodes_n_nodes(:)
  integer, allocatable :: nodes_n(:,:)
  integer              :: m, j_nod, k_nod
  real, allocatable    :: z_new(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  call Profiler % Start('Insert_Buildings')

  !--------------------------------!
  !   A couple of initial checks   !
  !--------------------------------!
  Assert(allocated(Grid % old_c))
  Assert(allocated(Grid % new_c))

  has_top = .false.
  do bc = Boundary_Regions()
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name(1:3) .eq. 'TOP') has_top = .true.
  end do
  if(.not. has_top) then
    call Message % Error(72,                                                &
                'The grid does not have a boundary condition called   ' //  &
                '"top".  Withouot it, you cannot define a computation ' //  &
                'domain for simulation of city flows \n \n            ' //  &
                'This error is critical.  Exiting!',                        &
                file=__FILE__, line=__LINE__)
  end if

  !----------------------------------------------------------!
  !                                                          !
  !   Phase I: Read the ground definition from an STL file   !
  !                                                          !
  !----------------------------------------------------------!

  print *, '#=============================================================' // &
           '============'
  print *, '# Enter the name of the ASCII STL ground file (with extenstion):'
  print *, '#-------------------------------------------------------------' // &
           '------------'
  name_in = File % Single_Word_From_Keyboard()
  call File % Open_For_Reading_Ascii(name_in, fu)

  !--------------------------!
  !   Count all the facets   !
  !--------------------------!
  rewind(fu)
  f = 0       ! initialize facet counter
  do
    call File % Read_Line(fu)
    if(Line % tokens(1) .eq. 'endsolid') exit
    if(Line % tokens(1) .eq. 'facet') f = f + 1
  end do
  print '(a38,i9)', '# Number of facets on the ground:    ', f

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  n_ground_facets = f
  allocate(ground_x(n_ground_facets, 3))
  allocate(ground_y(n_ground_facets, 3))
  allocate(ground_z(n_ground_facets, 3))

  !---------------------------------!
  !   Read all vertex coordinates   !
  !---------------------------------!
  area = 0.0
  rewind(fu)
  f = 0       ! re-initialize facet counter
  do
    call File % Read_Line(fu)
    if(Line % tokens(1) .eq. 'endsolid') exit
    if(Line % tokens(1) .eq. 'facet') then
      f = f + 1
      call File % Read_Line(fu)    ! 'outer loop'
      do i_ver = 1, 3
        call File % Read_Line(fu)  ! 'vertex 1, 2 and 3'
        read(Line % tokens(2), *) ground_x(f, i_ver)
        read(Line % tokens(3), *) ground_y(f, i_ver)
        read(Line % tokens(4), *) ground_z(f, i_ver)
      end do
      call File % Read_Line(fu)    ! 'endloop'
      area = area + Math % Triangle_Area(ground_x(f,1),  ground_y(f,1),  &
                                         ground_x(f,2),  ground_y(f,2),  &
                                         ground_x(f,3),  ground_y(f,3))
    end if
  end do
  print *, '# Read all ground facets!'

  x_max = maxval(maxval(ground_x,1));  x_min = minval(minval(ground_x,1))
  y_max = maxval(maxval(ground_y,1));  y_min = minval(minval(ground_y,1))
  z_max = maxval(maxval(ground_z,1));  z_min = minval(minval(ground_z,1))

  print '(a38,6f12.4)', '# Bounding box: (corner 1)           ',  &
                                                       x_min, y_min, z_min
  print '(a38,6f12.4)', '#               (corner 2)           ',  &
                                                       x_max, y_max, z_max
  print '(a38, f12.4)', '# Ground area from bounding box:     ',  &
                                          (x_max - x_min) * (y_max - y_min)
  print '(a38, f12.4)', '# Ground area as a sum of facets:    ', area
  Assert(Math % Approx_Real((x_max-x_min)*(y_max-y_min), area))

  !----------------------------------------------------------!
  !                                                          !
  !   Phase II: Sort nodes and cells, align where possible   !
  !                                                          !
  !----------------------------------------------------------!

  ! Sort nodes in columns in z direction
  call Grid % Sort_Nodes_By_Coordinates()

  ! Search for nodal clusters
  call Grid % Search_Coordinate_Clusters(nodal = .true.,   &
                               enforce_uniform = .false.,  &
                                            nz = a_node_layers)
  print '(a38,i9)', '# Absolute number of node layers:    ', a_node_layers
  Assert(mod(Grid % n_nodes, a_node_layers) .eq. 0)

  ! Estimate number of horizontal layers again (to check if sorting worked)
  do n = 1, Grid % n_nodes
    if(Grid % zn(n+1) < Grid % zn(n)) then
      n_node_layers = n
      exit
    end if
  end do
  print '(a38,i9)', '# Naive number of node layers:       ', n_node_layers

  ! Align node coordinates
  do n = 1, Grid % n_nodes, n_node_layers  ! these are ground nodes
    do nu = n + 1, n + n_node_layers - 1
      Grid % xn(nu) = Grid % xn(n)
      Grid % yn(nu) = Grid % yn(n)
    end do
  end do

  ! Search for celular clusters
  call Grid % Calculate_Cell_Centers()
  call Grid % Sort_Cells_By_Coordinates()
  call Grid % Search_Coordinate_Clusters(nodal = .false.,  &
                               enforce_uniform = .false.,  &
                                            nz = a_cell_layers)
  print '(a38,i9)', '# Absolute number of cell layers:    ', a_cell_layers
  Assert(mod(Grid % n_cells, a_cell_layers) .eq. 0)

  ! Estimate number of horizontal layers again (to check if sorting worked)
  do c = 1, Grid % n_cells
    if(Grid % zc(c+1) < Grid % zc(c)) then
      n_cell_layers = c
      exit
    end if
  end do
  print '(a38,i9)', '# Naive number of cell layers:       ', n_cell_layers

  ! Cells have been sorted before calling this function, which should hold
  Assert(a_cell_layers .eq. n_cell_layers)

  ! These browses only through cells "on the ground"
  do c = 1, Grid % n_cells, n_cell_layers
    Assert(any(Grid % cells_bnd_region(1:6, c) .ne. 0))
  end do

  !------------------------------------------------------------------!
  !                                                                  !
  !   Phase III: Find cells in buildingsi, chimneys and porosities   !
  !                                                                  !
  !------------------------------------------------------------------!

  !--------------------------------!
  !   Treat them all in one loop   !
  !--------------------------------!
  n_buildings  = 0
  n_porosities = 0
  allocate(cell_in_building(Grid % n_cells));  cell_in_building(:) = 0
  allocate(cell_in_chimney (Grid % n_cells));  cell_in_chimney (:) = 0
  Assert(allocated(Grid % por));               Grid % por      (:) = 0

  do run = FOR_BUILDINGS, FOR_POROSITIES
  ! Browse through cells "on the ground"
  do c = 1, Grid % n_cells, n_cell_layers
    z_min = HUGE
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      z_min = min(z_min, Grid % zn(n))
    end do
    do dir = 1, 6
      if(Grid % cells_bnd_region(dir, c) .ne. 0) then
        bc_name = trim(Grid % region % name(Grid % cells_bnd_region(dir, c)))
        call String % To_Upper_Case(bc_name)
        if(run .eq. FOR_BUILDINGS  .and. bc_name(1:8) .eq. 'BUILDING' .or.  &
           run .eq. FOR_CHIMNEYS   .and. bc_name(1:7) .eq. 'CHIMNEY'  .or.  &
           run .eq. FOR_POROSITIES .and. bc_name(1:8) .eq. 'POROSITY') then
          if(run .eq. FOR_BUILDINGS .or. run .eq. FOR_POROSITIES) then
            if(len(trim(bc_name)) .lt. 14) then
              call Message % Error(80,                                     &
                 'This version of Convert % Insert_Building, expects ' //  &
                 'the building number as the part of the name.  You  ' //  &
                 'don''t seem to have provided it.  Building''s name ' //  &
                 'is building_HHH (where HHH is three digits for the ' //  &
                 'height, but the program expects building_HHH_RRR   ' //  &
                 'where RRR is the rank (number) of the building.    ' //  &
                 ' \n \n Please fix the .geo script,                 ' //  &
                 're-run the GMSH, and try again.                    ',    &
                 file=__FILE__, line=__LINE__)
            end if
            read(bc_name(10:12), *) height
            read(bc_name(14:16), *) region  ! building's or porosity's rank
          else if(run .eq. FOR_CHIMNEYS) then
            read(bc_name( 9:11), *) height
            read(bc_name(13:15), *) region  ! building on which it is mounted
          end if
          do cu = c, c + n_cell_layers - 1
            z_max = -HUGE
            do i_nod = 1, abs(Grid % cells_n_nodes(cu))
              n = Grid % cells_n(i_nod, cu)
              z_max = max(z_max, Grid % zn(n))
            end do
            if(z_max - z_min > height + MICRO) goto 1

            if(run .eq. FOR_BUILDINGS) then
              cell_in_building(cu) = region  ! mark cell as in building
              n_buildings = max(n_buildings, region)
            end if

            if(run .eq. FOR_CHIMNEYS) then
              cell_in_chimney(cu) = region   ! mark cell as in this chimney
            end if

            if(run .eq. FOR_POROSITIES) then
              Grid % por(cu) = region        ! mark cell in porosity
              n_porosities = max(n_porosities, region)
            end if

          end do
        end if
      end if
    end do
1   continue  ! next near-wall cell
  end do
  end do

  ! Count the number of cells in buildings (this is for checking only)
  cnt = 0
  do c = 1, Grid % n_cells
    if(cell_in_building(c) .gt. 0) cnt = cnt + 1
  end do
  print '(a38,i9)', '# Number of buildings:               ', n_buildings
  print '(a38,i9)', '# Number of cells in buildings:      ', cnt

  ! Count the number of cells in porosities (this is for checking only)
  cnt = 0
  do c = 1, Grid % n_cells
    if(Grid % por(c) > 0) cnt = cnt + 1
  end do
  print '(a38,i9)', '# Number of porosities:              ', n_porosities
  print '(a38,i9)', '# Number of cells in porosities:     ', cnt

  !-----------------------------!
  !   Mark nodes on buildings   !
  !-----------------------------!
  allocate(node_on_building(Grid % n_nodes))
  node_on_building(:) = 0
  do c = 1, Grid % n_cells
    if(cell_in_building(c) .gt. 0) then
      do i_nod = 1, Grid % cells_n_nodes(c)  ! mark node as on building
        n = Grid % cells_n(i_nod, c)
        node_on_building(n) = cell_in_building(c)
      end do
    end if
  end do

  !----------------------------!
  !   Mark nodes on chimneys   !
  !----------------------------!
  allocate(node_on_chimney(Grid % n_nodes))
  node_on_chimney(:) = 0
  do c = 1, Grid % n_cells
    if(cell_in_chimney(c) > 0) then
      do i_nod = 1, Grid % cells_n_nodes(c)  ! mark node as on chimney
        n = Grid % cells_n(i_nod, c)
        node_on_chimney(n) = cell_in_chimney(c)
      end do
    end if
  end do

  !---------------------------------------------------------------!
  !                                                               !
  !   Phase IV: Place grid on the terrain defind in an STL file   !
  !                                                               !
  !---------------------------------------------------------------!

  print *, '#=================================================================='
  print *, '# Placing grid on the ground.  This should be much faster now.     '
  print *, '#------------------------------------------------------------------'

  ! Find maximum z (sky coordinate)
  z_min =  HUGE
  z_max = -HUGE
  do n = 1, Grid % n_nodes
    z_min = min(Grid % zn(n), z_min)
    z_max = max(Grid % zn(n), z_max)
  end do
  print '(a38,2f12.4)', '# Minimum and maximum z coordinate:  ', z_min, z_max

  !----------------------------------------------!
  !   For each node, find corresponding ground   !
  !     facet, and approximate z coordinate      !
  !----------------------------------------------!
  do n = 1, Grid % n_nodes, n_node_layers  ! this will only go ...
                                           ! ... through ground nodes
    do f = 1, n_ground_facets

      area = Math % Triangle_Area(ground_x(f,1),  ground_y(f,1),   &
                                  ground_x(f,2),  ground_y(f,2),   &
                                  ground_x(f,3),  ground_y(f,3))

      area_1 = Math % Triangle_Area(Grid % xn(n),  Grid % yn(n),   &
                                    ground_x(f,2), ground_y(f,2),  &
                                    ground_x(f,3), ground_y(f,3))

      area_2 = Math % Triangle_Area(Grid % xn(n),  Grid % yn(n),   &
                                    ground_x(f,3), ground_y(f,3),  &
                                    ground_x(f,1), ground_y(f,1))

      area_3 = Math % Triangle_Area(Grid % xn(n),  Grid % yn(n),   &
                                    ground_x(f,1), ground_y(f,1),  &
                                    ground_x(f,2), ground_y(f,2))

      if((area_1+area_2+area_3 - NANO) < area) then
        z_ground = (area_1 / area) * ground_z(f,1)  &
                 + (area_2 / area) * ground_z(f,2)  &
                 + (area_3 / area) * ground_z(f,3)
        do nu = n, n + n_node_layers - 1
          z_from_top    = z_max - Grid % zn(nu)
          z_scale       = (z_max - z_ground) / (z_max - z_min)
          Grid % zn(nu) = z_max - z_from_top * z_scale
        end do
      end if
    end do
  end do

  !------------------------------------------------!
  !                                                !
  !   Phase V: Flatten the tops of the buildings   !
  !                                                !
  !------------------------------------------------!

  ! Count number of nodes in x-y plane (orthogonal projection)
  cnt = 0
  do n = 1, Grid % n_nodes, n_node_layers  ! this will only go ...
    cnt = cnt + 1
  end do

  ! With this count, allocate memory for vertical (z) lines
  ! (This is clearly an overkill, but can optimize later)
  allocate(buildings_n_ground_nodes(n_buildings))
  allocate(buildings_ground_nodes  (n_buildings, cnt))
  buildings_n_ground_nodes(:) = 0
  buildings_ground_nodes(:,:) = 0

  ! Find all the vertical (z) lines in the buildings
  do n = 1, Grid % n_nodes, n_node_layers  ! this will only go ...
    b = node_on_building(n)
    if(b .gt. 0) then
      buildings_n_ground_nodes(b) = buildings_n_ground_nodes(b) + 1
      buildings_ground_nodes(b, buildings_n_ground_nodes(b)) = n
    end if
  end do

  !-----------------------------------------!
  !   Browse through all the buildings to   !
  !   flatten the lines from the roof up    !
  !-----------------------------------------!
  do b = 1, n_buildings
    do lay = 1, n_node_layers - 2  ! 0 would be ground, don't touch it
      z_avg = 0.0
      do i_nod = 1, buildings_n_ground_nodes(b)
        n = buildings_ground_nodes(b, i_nod)     ! ground node
        z_avg = z_avg + Grid % zn(n+lay)         ! z at lay
      end do
      z_avg = z_avg / buildings_n_ground_nodes(b)
      do i_nod = 1, buildings_n_ground_nodes(b)
        n = buildings_ground_nodes(b, i_nod)     ! ground node
        Grid % zn(n+lay) = z_avg                 ! z at lay
      end do
    end do
  end do

  !---------------------------------------------------!
  !   Smooth out the lines along the building a bit   !
  !---------------------------------------------------!
  print '(a)', ' # Smoothing out the grid lines along the buildings'
  do iter = 1, ITERS
    do b = 1, n_buildings

      ! Search for the top (layer) of this bulding)
      do lay = 1, n_node_layers - 1  ! 0 would be ground, don't touch it
        i_nod = 1
        n = buildings_ground_nodes(b, i_nod)     ! ground node
        Assert(node_on_building(n) .eq. b)
        if(node_on_building(n+lay+lay) .eq. 0) then
          top = lay
          if(iter .eq. 1) then
            print '(a,i4,a,i4)', ' # Top layer on bulding', b, ' is: ', top
            Assert(node_on_building(n+top) .eq. b)
          end if
          exit
        end if
      end do    ! layer

      ! Smooth the vertical lines on the buildings
      do lay = 1, top - 1
        do i_nod = 1, buildings_n_ground_nodes(b)
          n = buildings_ground_nodes(b, i_nod)     ! ground node
          Grid % zn(n+lay) = 0.5 * (  Grid % zn(n+lay-1)  &
                                    + Grid % zn(n+lay+1))
          Assert(node_on_building(n+lay+1) .eq. b)
        end do  ! vertical line in the building
      end do    ! layer
    end do      ! building
  end do        ! iteration

  !-----------------------------------------------------------!
  !                                                           !
  !   Phase VI: Manage boundary condition names and numbers   !
  !                                                           !
  !-----------------------------------------------------------!

  ! Find ground b.c. number (used to eliminate BUILDING_000 and POROSITY_...)
  do bc = Boundary_Regions()
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name .eq. 'GROUND') then
      ground_bc = bc
    end if
  end do

  ! Eliminate building_000 and porosity_... b.c.
  ! (Where on earth could building_000 come from?)
  do bc = Boundary_Regions()
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name(1:12) .eq. 'BUILDING_000' .or.  &
       bc_name(1:8)  .eq. 'POROSITY') then
      do c = 1, Grid % n_cells
        do dir = 1, 6
          if(Grid % cells_bnd_region(dir, c) .eq. bc) then
            Grid % cells_bnd_region(dir, c) = ground_bc
          end if
        end do
      end do
    end if
  end do

  ! Counts non-building and non-porosity boundary conditions
  cnt   = 0
  do bc = Boundary_Regions()
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name(1:8) .ne. 'BUILDING' .and.  &
       bc_name(1:8) .ne. 'POROSITY') then
      cnt = cnt + 1
      Grid % region % name(cnt) = Grid % region % name(bc)
      do c = 1, Grid % n_cells
        do dir = 1, 6
          if(Grid % cells_bnd_region(dir, c) .eq. bc) then
            Grid % cells_bnd_region(dir, c) = cnt
          end if
        end do
      end do
    end if
  end do
  Grid % n_bnd_regions = cnt

  ! Shift all by one up (to be able to insert building walls as first)
  if(n_buildings .gt. 0) then
    do bc = Grid % n_bnd_regions, 1, -1
      Grid % region % name(bc+1) = Grid % region % name(bc)
      do c = 1, Grid % n_cells
        do dir = 1, 6
          if(Grid % cells_bnd_region(dir, c) .eq. bc) then
            Grid % cells_bnd_region(dir, c) = bc+1
          end if
        end do
      end do
    end do
    Grid % n_bnd_regions = Grid % n_bnd_regions + 1

    ! Add first boundary condition for walls
    Grid % region % name(1) = 'BUILDING_WALLS'
  end if

  !------------------------------------------------------------!
  !                                                            !
  !   Phase VII: Store only cells which are not in buildings   !
  !                                                            !
  !------------------------------------------------------------!

  !-----------------------------------------------------!
  !   Store cells' nodes and boundary conditons again   !
  !-----------------------------------------------------!
  allocate(i_work_1(   Grid % n_cells));  i_work_1(:)   = 0
  allocate(i_work_2(8, Grid % n_cells));  i_work_2(:,:) = 0
  allocate(i_work_3(6, Grid % n_cells));  i_work_3(:,:) = 0
  allocate(i_work_4(   Grid % n_cells));  i_work_4(:)   = 0

  do c = 1, Grid % n_cells
    i_work_1(c)      = Grid % cells_n_nodes(c)
    i_work_2(1:8, c) = Grid % cells_n(1:8, c)
    i_work_3(1:6, c) = Grid % cells_bnd_region(1:6, c)
    i_work_4(c)      = Grid % por(c)
  end do

  !------------------------------------------------------!
  !   Renumber the cells without the cells in building   !
  !------------------------------------------------------!
  Grid % new_c(:) = 0
  cnt             = 0
  do c = 1, Grid % n_cells
    if(cell_in_building(c) == 0 .and. cell_in_chimney(c) == 0) then
      cnt = cnt + 1
      Grid % new_c(c) = cnt
    end if
  end do
  print '(a38,i9)', '# Old number of cells:               ', Grid % n_cells
  print '(a38,i9)', '# Number of cells without buildings: ', cnt
  print '(a38,i9)', '# Number of excluded cells:          ', Grid % n_cells-cnt

  !--------------------------------------!
  !   Information on cell connectivity   !
  !--------------------------------------!
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then
      Grid % cells_n_nodes        (Grid % new_c(c)) = i_work_1(c)
      Grid % cells_n         (1:8, Grid % new_c(c)) = i_work_2(1:8, c)
      Grid % cells_bnd_region(1:6, Grid % new_c(c)) = i_work_3(1:6, c)
      Grid % por                  (Grid % new_c(c)) = i_work_4(c)
    end if
  end do
  Grid % n_cells = cnt

  !-------------------------------------------------------------------------!
  !                                                                         !
  !   Phase VIII: Insert new boundary faces on building and chimney walls   !
  !                                                                         !
  !-------------------------------------------------------------------------!
  if(n_buildings .gt. 0) then

    do c = 1, Grid % n_cells

      if(Grid % cells_n_nodes(c) .eq. 4) fn = TET
      if(Grid % cells_n_nodes(c) .eq. 5) fn = PYR
      if(Grid % cells_n_nodes(c) .eq. 6) fn = WED
      if(Grid % cells_n_nodes(c) .eq. 8) fn = HEX

      do dir = 1, 6
        if(Grid % cells_bnd_region(dir, c) .eq. 0) then

          n_f_nod    = 0
          f_nod(1:4) = -1

          ! Fill up face nodes
          do i_nod = 1, 4
            if(fn(dir, i_nod) > 0) then
              f_nod(i_nod) = Grid % cells_n(fn(dir, i_nod), c)
              n_f_nod   = n_f_nod + 1
            end if
          end do

          if( n_f_nod > 0 ) then

            ! Qadrilateral face
            if(f_nod(4) > 0) then
              if( (node_on_building(f_nod(1)) .gt. 0) .and.  &
                  (node_on_building(f_nod(2)) .gt. 0) .and.  &
                  (node_on_building(f_nod(3)) .gt. 0) .and.  &
                  (node_on_building(f_nod(4)) .gt. 0) ) then
                Grid % cells_bnd_region(dir, c) = 1
              end if
            else
              if( (node_on_building(f_nod(1)) .gt. 0) .and.  &
                  (node_on_building(f_nod(2)) .gt. 0) .and.  &
                  (node_on_building(f_nod(3)) .gt. 0) ) then
                Grid % cells_bnd_region(dir, c) = 1
              end if
            end if
          end if
        end if
      end do
    end do

  end if

  do bc = Boundary_Regions()
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name(1:7) .eq. 'CHIMNEY') then

      read(bc_name(13:15), *) region  ! this is, in essence, building number
                                      ! (4 years l8r: is it chimney number?)
      do c = 1, Grid % n_cells

        if(Grid % cells_n_nodes(c) .eq. 4) fn = TET
        if(Grid % cells_n_nodes(c) .eq. 5) fn = PYR
        if(Grid % cells_n_nodes(c) .eq. 6) fn = WED
        if(Grid % cells_n_nodes(c) .eq. 8) fn = HEX

        do dir = 1, 6
          if(Grid % cells_bnd_region(dir, c) .eq. 0) then

            n_f_nod    = 0
            f_nod(1:4) = -1

            ! Fill up face nodes
            do i_nod = 1, 4
              if(fn(dir, i_nod) > 0) then
                f_nod(i_nod) = Grid % cells_n(fn(dir, i_nod), c)
                n_f_nod   = n_f_nod + 1
              end if
            end do

            if( n_f_nod > 0 ) then

              ! Qadrilateral face
              if(f_nod(4) > 0) then
                if( node_on_chimney(f_nod(1)) .eq. region .and.  &
                    node_on_chimney(f_nod(2)) .eq. region .and.  &
                    node_on_chimney(f_nod(3)) .eq. region .and.  &
                    node_on_chimney(f_nod(4)) .eq. region ) then
                  Grid % cells_bnd_region(dir, c) = bc
                end if
              ! Triangular face
              else
                if( node_on_chimney(f_nod(1)) .eq. region .and.  &
                    node_on_chimney(f_nod(2)) .eq. region .and.  &
                    node_on_chimney(f_nod(3)) .eq. region ) then
                  Grid % cells_bnd_region(dir, c) = bc
                end if
              end if
            end if
          end if
        end do
      end do

    end if

  end do  ! chimney

  !----------------------------------------------------!
  !                                                    !
  !   Phase IX: Smooth the grid lines in z direction   !
  !                                                    !
  !----------------------------------------------------!

  print '(a)', ' # Smoothing the grid out a little bit'

  ! This is a wee-bit memory intensive, better be on the safe side
  deallocate(cell_in_building)
  deallocate(cell_in_chimney)
  deallocate(node_on_building)
  deallocate(node_on_chimney)
  deallocate(buildings_n_ground_nodes)
  deallocate(buildings_ground_nodes)
  deallocate(i_work_1)
  deallocate(i_work_2)
  deallocate(i_work_3)
  deallocate(i_work_4)

  allocate(nodes_n_nodes(Grid % n_nodes));  nodes_n_nodes(:) = 0
  allocate(nodes_n(MAX_N,Grid % n_nodes));  nodes_n    (:,:) = 0
  allocate(z_new   (Grid % n_nodes))

  !-----------------------!
  !   Connect the nodes   !
  !-----------------------!
  do c = 1, Grid % n_cells                        ! through cells
    do i_nod = 1, abs(Grid % cells_n_nodes(c))    ! nodes of the cell
      n = Grid % cells_n(i_nod, c)                ! first node
      do j_nod = 1, abs(Grid % cells_n_nodes(c))  ! nodes of the cell
        if(j_nod .ne. i_nod) then
          m = Grid % cells_n(j_nod, c)            ! second node
          do k_nod = 1, nodes_n_nodes(n)
            if(nodes_n(k_nod, n) .eq. m) goto 4
          end do
          nodes_n_nodes(n) = nodes_n_nodes(n) + 1
          Assert(nodes_n_nodes(n) .le. MAX_N)
          nodes_n(nodes_n_nodes(n), n) = m
4         continue
        end if
      end do
    end do
  end do

  ! Sort them, for what it's worth
  do n = 1, Grid % n_nodes
    if(nodes_n_nodes(n) .gt. 0) then  ! avoid nodes in cells which have ...
                                      ! ... been deleted (it wasn't checked)
      call Sort % Int_Array(nodes_n(1:nodes_n_nodes(n), n))
    end if
  end do

  !-----------------------------------------!
  !   Exclude nodes which should not move   !
  !-----------------------------------------!
  do c = 1, Grid % n_cells
    if(Grid % cells_n_nodes(c) .eq. 4) fn = TET
    if(Grid % cells_n_nodes(c) .eq. 5) fn = PYR
    if(Grid % cells_n_nodes(c) .eq. 6) fn = WED
    if(Grid % cells_n_nodes(c) .eq. 8) fn = HEX
    do dir = 1, 6
      if(Grid % cells_bnd_region(dir, c) .ne. 0) then
        bc_name = trim(Grid % region % name(Grid % cells_bnd_region(dir, c)))
        call String % To_Upper_Case(bc_name)
        if(bc_name(1:8) .eq. 'BUILDING' .or.  &
           bc_name(1:7) .eq. 'CHIMNEY'  .or.  &
           bc_name(1:6) .eq. 'GROUND'   .or.  &
           bc_name(1:3) .eq. 'TOP') then
          do i_nod = 1, 4                ! faces' nodes
            if(fn(dir, i_nod) > 0) then  ! this node exists
              n = Grid % cells_n(fn(dir, i_nod), c)
              nodes_n_nodes(n) = 0       ! brute but practical way to exclude it
            end if
          end do
        end if
      end if
    end do
  end do

  !-------------------------------------!
  !   Smooth in the z direction a bit   !
  !-------------------------------------!
  do iter = 1, ITERS

    ! Compute averaged z coordinate for each node
    do n = 1, Grid % n_nodes
      if(nodes_n_nodes(n) > 0) then
        z_avg = 0.0
        do i_nod = 1, nodes_n_nodes(n)
          z_avg = z_avg + Grid % zn(nodes_n(i_nod, n))
        end do
        z_avg = z_avg / real(nodes_n_nodes(n))
        z_new(n) = z_avg
      end if
    end do  ! through nodes

    ! Perform the actual smoothing
    do n = 1, Grid % n_nodes
      if(nodes_n_nodes(n) > 0) then

        ! No matter how small, some under-relaxation is needed here
        Grid % zn(n) = 0.1 * Grid % zn(n) + 0.9 * z_new(n)
      end if
    end do  ! through nodes
  end do    ! iterations

  call Profiler % Stop('Insert_Buildings')

  end subroutine
