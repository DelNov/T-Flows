!==============================================================================!
  subroutine Insert_Buildings(Convert, Grid)
!------------------------------------------------------------------------------!
!>  This soubroutine was developed to import mesh with buildings, and place
!>  them both on a ground defined with an STL file.  This subroutine was
!>  traditionally not a part of the standard T-Flows distribution, but after
!>  publishing a couple of papers with it, we decided to include it.  It is
!>  rather long and often chaning, and it will probably stay like that in the
!>  foreseeable future.
!>
!>  The vertical distribution of the nodes is taken from the (flat) grid, so
!>  it is defined in the .geo script.  A few more parameters can be given in
!>  an optional file <problem_name>.ins (e.g. urban.city.ins for the grid
!>  urban.city.msh), one key and value per line, which a .geo script can
!>  write with Printf("ROUGHNESS_ALLOWED %g", 0.5) > "urban.city.ins"; :
!>
!>  ROUGHNESS_ALLOWED   [m]  terrain nodes farther than that from the plane
!>                           through their neighbors are moved towards it
!>                           (default: terrain is kept as in the STL file)
!>  HEIGHT_ABOVE_ROOFS  [m]  zone above the tallest building in which the
!>                           deformations caused by roofs are absorbed
!>                           (default: max(0.5 H_max, 5 dz at the roof))
!>  Z_SMOOTH_ITERS           iterations of smoothing in z direction at the
!>                           end (default: environment variable
!>                           CITY_Z_SMOOTH_ITERS if set, 24 otherwise)
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

  ! Parameters from the optional file <problem_name>.ins
  character(len=SL) :: ins_name
  logical           :: ins_found, reached_end
  integer           :: n_smooth, env_stat
  character(len=16) :: env_val
  real              :: h_above_ins

  ! Variables defining the vertical distribution of the nodes: the zone
  ! holding the buildings and a buffer above the tallest of them keeps
  ! the layers of the flat grid, the layers above it are scaled to reach
  ! the top of the domain
  real              :: dz_bld       ! layer thickness at the tallest roof
  real              :: h_above      ! zone above the tallest building
  real              :: h_zone       ! height of the zone of buildings
  integer           :: k_zone       ! node level at the top of the zone
  integer           :: n_lay
  real              :: z_flat_bot, z_flat_top, max_obstacle_height
  real              :: height_bc, l_up, s_top, scale_up, dev_max
  real              :: r_min, r_max
  integer           :: ios
  real,    parameter :: SAFE_ROUGH = 0.75 ! safe first cell / max roughness
  integer, parameter :: MAX_NB     = 16   ! max neighbors of a ground node
  real                 :: h_first         ! height of the first cell
  real                 :: h_first_safe    ! its safe value for the terrain
  real                 :: dz
  real,    allocatable :: s_lay(:)        ! heights of levels above ground
  real,    allocatable :: col_zg(:)       ! elevation of ground in columns
  real,    allocatable :: col_rg(:)       ! roughness of ground in columns
  real,    allocatable :: col_pl(:)       ! plane through the neighbors
  real,    allocatable :: col_zg0(:)      ! ground elevation from the STL
  integer, allocatable :: col_srt(:)      ! smoothed columns, largest first
  real,    allocatable :: chg(:)          ! their changes of elevation
  integer, parameter   :: MAX_SMOOTH_ITERS = 100
  real                 :: rough_allow     ! largest roughness allowed
  integer, allocatable :: col_n_nb(:), col_nb(:,:)  ! neighbors on ground
  integer              :: col_rough
  real                 :: rough, rough_max, rough_avg, a_pl, det, dx, dy
  real                 :: s11, s12, s13, s22, s23, s33, r1, r2, r3

  ! New variables for phase III: searching for cells in buildings
  ! (It's more accurate doing this before placing grid on ground)
  integer, allocatable :: cell_in_building(:)
  integer, allocatable :: cell_in_chimney(:),  node_on_chimney(:)
  character(len=SL)    :: bc_name
  real                 :: height                  ! building height
  integer              :: run, region, dir, cu, cnt, bc
  integer              :: n_buildings
  integer              :: n_porosities            ! for porosities
  real, allocatable    :: bld_height(:)           ! nominal building heights

  ! New variables for phase IV: place grid on the ground defined in STL
  real, parameter :: MAX_SCALE = 2.0   ! warn if the layers above the zone
                                       ! of buildings are scaled more
  real    :: area_1, area_2, area_3, z_ground
  integer :: columns_unmapped
  logical :: ground_found

  ! New variables for phase V: flatten the tops of the buildings
  integer, parameter   :: MAX_PINS          = 4    ! max roofs in a column
  real,    parameter   :: MIN_WALL_FRACTION = 0.5  ! min uphill wall height,
                                                   ! relative to the nominal
  real,    parameter   :: MIN_DZ_RATIO      = 0.6  ! max squeezing of layers
  integer              :: b, lay, top, iter, col_top
  integer              :: col, n_cols, n_roof, n_fixed, i, j, lay_a, lay_b
  integer, allocatable :: column_seen(:)
  integer, allocatable :: col_n_pins(:), col_pin_lay(:,:)
  real,    allocatable :: col_pin_z(:,:), z_old(:)
  real                 :: z_avg, z_roof, zg_min, zg_max, zg_avg
  real                 :: z_a, z_b, dz_min, ratio, h_meshed
  real                 :: ratio_min, ratio_max
  logical              :: consistent_roof

  ! New variables for phase VI: manage boundary condition names and numbers
  integer :: ground_bc

  ! New variables for phase VII: insert new boundary faces
  integer :: fn(6,4), n_f_nod, f_nod(4), s_nod(4)
  integer :: c1, c2, s
  logical :: c1_on, c2_on

  ! New variables for phase VIII: store only cells which are not in buildings
  integer, allocatable :: i_work_1(:)
  integer, allocatable :: i_work_2(:,:)
  integer, allocatable :: i_work_3(:,:)
  integer, allocatable :: i_work_4(:)

  ! New variables for phase IX: smooth the grid lines in z direction
  integer, allocatable :: nodes_n_nodes(:)
  integer, allocatable :: nodes_n(:,:)
  integer              :: m, j_nod, k_nod
  real, allocatable    :: z_new(:)
  real, allocatable    :: z_ref(:)  ! node elevations before the smoothing
  real, parameter      :: MIN_SMOOTH_DZ = 0.5  ! max squeezing of the layers
                                               ! allowed during smoothing
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
    call Message % Error(72,                                               &
                'The grid does not have a boundary condition called  ' //  &
                '"top".  Without it, you cannot define a computation ' //  &
                'domain for simulation of city flows \n \n           ' //  &
                'This error is critical.  Exiting!',                       &
                file=__FILE__, line=__LINE__)
  end if

  !----------------------------------------------------------!
  !                                                          !
  !   Phase I: Read the ground definition from an STL file   !
  !                                                          !
  !----------------------------------------------------------!

  print *, '#============================================================' // &
           '============'
  print *, '# Enter the name of the ASCII STL ground file (with extension):'
  print *, '#------------------------------------------------------------' // &
           '------------'
  name_in = File % Single_Word_From_Keyboard("# STL ground file (with ext.)")
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
  ! The STL needn't be a rectangle (it can also be larger than the grid),
  ! columns of the grid not covered by it are checked in Phase II b
  if(.not. Math % Approx_Real((x_max-x_min)*(y_max-y_min), area)) then
    print '(a)', ' # WARNING: The STL ground is not a rectangle aligned with'
    print '(a)', ' # the coordinate axes (area of facets differs from the one'
    print '(a)', ' # of the bounding box).'
  end if

  !------------------------------------------------------!
  !   Parameters from the optional file <problem>.ins    !
  !------------------------------------------------------!
  rough_allow = HUGE   ! terrain kept as in the STL file
  h_above_ins = -HUGE  ! default computed below
  n_smooth    = ITERS
  call get_environment_variable('CITY_Z_SMOOTH_ITERS', env_val,  &
                                status=env_stat)
  if(env_stat .eq. 0) read(env_val, *) n_smooth

  ins_name = Grid % name
  call String % To_Lower_Case(ins_name)
  ins_name = trim(ins_name) // '.ins'
  inquire(file=trim(ins_name), exist=ins_found)
  if(ins_found) then
    print '(a)', ' # Reading parameters from: ' // trim(ins_name)
    call File % Open_For_Reading_Ascii(ins_name, fu)
    do
      call File % Read_Line(fu, reached_end=reached_end)
      if(reached_end) exit
      if(Line % n_tokens < 2) cycle
      call String % To_Upper_Case(Line % tokens(1))
      select case(Line % tokens(1))
        case('ROUGHNESS_ALLOWED')
          read(Line % tokens(2), *) rough_allow
        case('HEIGHT_ABOVE_ROOFS')
          read(Line % tokens(2), *) h_above_ins
        case('Z_SMOOTH_ITERS')
          read(Line % tokens(2), *) height_bc  ! %g can write 0 as 0.0
          n_smooth = nint(height_bc)
        case default
          print '(a)', ' # WARNING: unknown key in ' // trim(ins_name) //  &
                       ': ' // trim(Line % tokens(1))
      end select
    end do
    close(fu)
  end if

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

  !-------------------------------------------------------------!
  !                                                             !
  !   Phase II b: Define the vertical distribution of nodes     !
  !                                                             !
  !   The heights of the node levels above the ground are taken !
  !   from the flat grid (they are defined in the .geo script). !
  !   They are kept in the zone holding all the buildings and   !
  !   h_above over the tallest of them.  In that zone, roofs    !
  !   are flattened, and the deformation they cause is fully    !
  !   absorbed in it.  Above, the layers are scaled to reach    !
  !   the top of the domain.  Terrain roughness is evaluated    !
  !   (and optionally reduced) here as well.                    !
  !                                                             !
  !-------------------------------------------------------------!

  ! The lowest and highest levels of the ORIGINAL (flat) grid
  z_flat_bot = minval(Grid % zn(1:Grid % n_nodes))
  z_flat_top = maxval(Grid % zn(1:Grid % n_nodes))
  if(z_flat_top <= z_flat_bot) then
    call Message % Error(84, 'Invalid flat mesh: z_max must exceed z_min.', &
                         file=__FILE__, line=__LINE__)
  end if

  ! Height of the tallest obstacle, from the boundary condition names
  max_obstacle_height = 0.0
  do bc = Boundary_Regions()
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    ios = 1
    if(bc_name(1:9) .eq. 'BUILDING_' .or. &
       bc_name(1:9) .eq. 'POROSITY_') then
      read(bc_name(10:12), *, iostat=ios) height_bc
    else if(bc_name(1:8) .eq. 'CHIMNEY_') then
      read(bc_name(9:11), *, iostat=ios) height_bc
    end if
    if(ios .eq. 0) max_obstacle_height = max(max_obstacle_height, height_bc)
  end do
  print '(a38,f12.4)', '# Highest obstacle height [m]:        ', max_obstacle_height

  !----------------------------------------------!
  !   Elevation of the terrain in each column    !
  !----------------------------------------------!
  n_cols = Grid % n_nodes / n_node_layers
  allocate(col_zg(n_cols));  col_zg(:) = 0.0

  columns_unmapped = 0
  do col = 1, n_cols
    n = (col-1) * n_node_layers + 1  ! ground node
    ground_found = .false.
    do f = 1, n_ground_facets

      area = Math % Triangle_Area(ground_x(f,1), ground_y(f,1), &
                                  ground_x(f,2), ground_y(f,2), &
                                  ground_x(f,3), ground_y(f,3))
      if(area <= NANO) cycle  ! Degenerate ground triangle

      area_1 = Math % Triangle_Area(Grid % xn(n),  Grid % yn(n), &
                                    ground_x(f,2), ground_y(f,2), &
                                    ground_x(f,3), ground_y(f,3))
      area_2 = Math % Triangle_Area(Grid % xn(n),  Grid % yn(n), &
                                    ground_x(f,3), ground_y(f,3), &
                                    ground_x(f,1), ground_y(f,1))
      area_3 = Math % Triangle_Area(Grid % xn(n),  Grid % yn(n), &
                                    ground_x(f,1), ground_y(f,1), &
                                    ground_x(f,2), ground_y(f,2))

      ! A point exactly on two STL triangles must not be mapped twice.
      if((area_1 + area_2 + area_3 - NANO) < area) then
        col_zg(col) = (area_1 / area) * ground_z(f,1) &
                    + (area_2 / area) * ground_z(f,2) &
                    + (area_3 / area) * ground_z(f,3)
        ground_found = .true.
        exit
      end if
    end do
    if(.not. ground_found) then
      columns_unmapped = columns_unmapped + 1
      if(columns_unmapped <= 5) then
        print '(a,2f14.4)', ' # WARNING: STL does not cover x,y: ', &
                           Grid % xn(n), Grid % yn(n)
      end if
    end if
  end do
  if(columns_unmapped > 0) then
    call Message % Error(88,                                       &
      'Some mesh columns lie outside the STL ground. Extend the STL.', &
      file=__FILE__, line=__LINE__)
  end if

  !-------------------------------------------------------------------!
  !   Roughness of the terrain, as the grid sees it.                  !
  !                                                                   !
  !   For each ground node, a plane is fitted (least squares) through !
  !   its neighbors on the ground, and roughness is the distance of   !
  !   the node from that plane.  A uniform slope is harmless for the  !
  !   terrain-following layers, but a bump (or a dip) of one node     !
  !   bends the near-ground layers on the scale of a single cell.     !
  !   If the first cell is thinner than such a bump, the near-ground  !
  !   layer of cells gets folded around it, and the simulations were !
  !   seen to diverge there.  The safe height of the first cell is    !
  !   therefore taken as a multiple of the largest roughness.  This   !
  !   holds for least squares gradients of pressure, Gauss theorem   !
  !   for pressure gradients was stable with much thinner cells.     !
  !-------------------------------------------------------------------!
  allocate(col_n_nb(n_cols));          col_n_nb(:)   = 0
  allocate(col_nb(MAX_NB, n_cols));    col_nb  (:,:) = 0
  do c = 1, Grid % n_cells, n_cell_layers  ! cells on the ground
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      if(mod(n-1, n_node_layers) .ne. 0) cycle  ! not on the ground
      col = (n-1) / n_node_layers + 1
      do j = 1, abs(Grid % cells_n_nodes(c))
        m = Grid % cells_n(j, c)
        if(mod(m-1, n_node_layers) .ne. 0) cycle
        i = (m-1) / n_node_layers + 1
        if(i .eq. col) cycle
        if(any(col_nb(1:col_n_nb(col), col) .eq. i)) cycle
        if(col_n_nb(col) .lt. MAX_NB) then
          col_n_nb(col) = col_n_nb(col) + 1
          col_nb(col_n_nb(col), col) = i
        end if
      end do
    end do
  end do

  allocate(col_rg (n_cols));  col_rg (:) = 0.0
  allocate(col_pl (n_cols));  col_pl (:) = 0.0
  allocate(col_zg0(n_cols));  col_zg0(:) = col_zg(:)

  !--------------------------------------------------------------!
  !   Optional local smoothing of the terrain: ground nodes      !
  !   farther than rough_allow (ROUGHNESS_ALLOWED in .ins file)  !
  !   from the plane through their neighbors are moved towards   !
  !   it, just to that distance.  Since that changes the planes  !
  !   of the neighbors too, it is repeated until no node exceeds !
  !   the allowed roughness.                                     !
  !--------------------------------------------------------------!
  do iter = 0, MAX_SMOOTH_ITERS

    rough_max = 0.0
    rough_avg = 0.0
    col_rough = 0
    do col = 1, n_cols
      col_rg(col) = 0.0
      col_pl(col) = col_zg(col)
      if(col_n_nb(col) < 3) cycle
      n = (col-1) * n_node_layers + 1
      ! Normal equations of the plane z = a + b*dx + c*dy around the node
      s11 = 0.0; s12 = 0.0; s13 = 0.0; s22 = 0.0; s23 = 0.0; s33 = 0.0
      r1  = 0.0; r2  = 0.0; r3  = 0.0
      do j = 1, col_n_nb(col)
        i  = col_nb(j, col)
        m  = (i-1) * n_node_layers + 1
        dx = Grid % xn(m) - Grid % xn(n)
        dy = Grid % yn(m) - Grid % yn(n)
        s11 = s11 + 1.0;    s12 = s12 + dx;     s13 = s13 + dy
        s22 = s22 + dx*dx;  s23 = s23 + dx*dy;  s33 = s33 + dy*dy
        r1  = r1  + col_zg(i)
        r2  = r2  + col_zg(i) * dx
        r3  = r3  + col_zg(i) * dy
      end do
      det = s11 * (s22*s33 - s23*s23)  &
          - s12 * (s12*s33 - s23*s13)  &
          + s13 * (s12*s23 - s22*s13)
      if(abs(det) < MICRO) cycle  ! neighbors on a line
      ! Only the value of the plane at the node (a) is needed
      a_pl = ( r1  * (s22*s33 - s23*s23)  &
             - s12 * (r2 *s33 - s23*r3 )  &
             + s13 * (r2 *s23 - s22*r3 ) ) / det
      rough = abs(col_zg(col) - a_pl)
      col_pl(col) = a_pl
      col_rg(col) = rough
      rough_avg = rough_avg + rough
      if(rough > rough_max) then
        rough_max = rough
        col_rough = col
      end if
    end do
    rough_avg = rough_avg / real(n_cols)

    ! Report the terrain as read from the STL file
    if(iter .eq. 0) then
      print '(a38,f12.4)', '# Average roughness of terrain [m]:   ', rough_avg
      print '(a38,f12.4)', '# Largest roughness of terrain [m]:   ', rough_max
      if(col_rough > 0) then
        n = (col_rough-1) * n_node_layers + 1
        print '(a38,3f12.2)', '# ... found at (x, y, z):             ',  &
                              Grid % xn(n), Grid % yn(n), col_zg(col_rough)
      end if
      do i = 1, 5
        dz = 0.25 * real(i) * rough_max
        print '(a,f7.3,a,i9)', ' # Ground nodes with roughness above ', dz,  &
                               ' m:', count(col_rg(:) > dz)
      end do
      print '(a38,f12.4)', '# Safe first cell for STL terrain [m]:',  &
                           SAFE_ROUGH * rough_max
      print '(a)', ' # (with least squares gradients of pressure; with Gauss'
      print '(a)', ' #  theorem for pressure gradients it can be thinner)'
      if(rough_allow .ge. HUGE) exit  ! no smoothing of the terrain
      print '(a38,f12.4)', '# Largest roughness allowed [m]:      ', rough_allow
    end if

    ! Move the nodes which are too rough (all at once, from the planes
    ! computed above, so that the result doesn't depend on numbering)
    cnt = 0
    do col = 1, n_cols
      if(col_rg(col) > rough_allow + NANO) then
        col_zg(col) = col_pl(col)  &
                    + sign(rough_allow, col_zg(col) - col_pl(col))
        cnt = cnt + 1
      end if
    end do
    if(cnt .eq. 0) exit
  end do

  !------------------------------------------!
  !   Report where and how much the terrain  !
  !   was smoothed, on screen and in a file  !
  !------------------------------------------!
  if(rough_allow < HUGE) then
    if(iter > MAX_SMOOTH_ITERS) then
      print '(a)', ' # WARNING: Smoothing of the terrain did not converge;'
      print '(a,f9.4,a)', ' # some nodes remain rougher than ',  &
                          rough_allow, ' m'
    end if

    ! Indices of the smoothed nodes; only the 20 largest changes are
    ! printed, so they are picked one by one instead of sorting them all
    cnt = count(abs(col_zg(:) - col_zg0(:)) > NANO)
    allocate(col_srt(cnt))
    allocate(chg(cnt))
    cnt = 0
    do col = 1, n_cols
      if(abs(col_zg(col) - col_zg0(col)) > NANO) then
        cnt = cnt + 1
        col_srt(cnt) = col
        chg(cnt) = abs(col_zg(col) - col_zg0(col))
      end if
    end do
    do i = 1, min(cnt, 20)
      j = i - 1 + maxloc(chg(i:cnt), 1)
      col = col_srt(j);  col_srt(j) = col_srt(i);  col_srt(i) = col
      dz  = chg(j);      chg(j)     = chg(i);      chg(i)     = dz
    end do

    print '(a)', ' #-----------------------------------------------------'
    print '(a,f9.4,a,i4,a)', ' # Terrain smoothed to roughness ', rough_allow,  &
                             ' m in ', iter, ' iterations'
    print '(a38,i9)', '# Ground nodes moved:                ', cnt
    if(cnt > 0) then
      print '(a38,f12.4)', '# Largest change of elevation [m]:   ',  &
                           abs(col_zg(col_srt(1)) - col_zg0(col_srt(1)))
      print '(a38,f12.4)', '# Average change of moved nodes [m]: ',  &
                           sum(abs(col_zg(col_srt(:)) - col_zg0(col_srt(:))))  &
                           / real(cnt)
      print '(a)', ' #          x           y     z (STL)   z (smooth)' //  &
                   '     change'
      do i = 1, min(cnt, 20)
        col = col_srt(i)
        n   = (col-1) * n_node_layers + 1
        print '(a,2f12.2,2f12.3,f11.3)', ' # ', Grid % xn(n), Grid % yn(n),  &
              col_zg0(col), col_zg(col), col_zg(col) - col_zg0(col)
      end do
      if(cnt > 20) print '(a,i9,a)', ' # ... and ', cnt - 20,  &
                                     ' smaller changes, see the file'

      ! All of them in a file (ParaView reads it with "Table To Points")
      call File % Open_For_Writing_Ascii('terrain_smoothing.csv', fu)
      write(fu, '(a)') 'x,y,z_stl,z_smooth,change'
      do i = 1, cnt
        col = col_srt(i)
        n   = (col-1) * n_node_layers + 1
        write(fu, '(2(f0.4,","),2(f0.4,","),f0.4)') Grid % xn(n),   &
              Grid % yn(n), col_zg0(col), col_zg(col),              &
              col_zg(col) - col_zg0(col)
      end do
      close(fu)
    end if

    ! Roughness after smoothing
    print '(a38,f12.4)', '# Largest roughness after smoothing: ', rough_max
    if(col_rough > 0) then
      n = (col_rough-1) * n_node_layers + 1
      print '(a38,3f12.2)', '# ... found at (x, y, z):             ',  &
                            Grid % xn(n), Grid % yn(n), col_zg(col_rough)
    end if
  end if

  h_first_safe = SAFE_ROUGH * rough_max

  !--------------------------------------------------------!
  !   Heights of the node levels above the ground, taken   !
  !   from the flat grid (all its columns must be equal)   !
  !--------------------------------------------------------!
  allocate(s_lay(0:n_node_layers-1))
  s_lay(0:n_node_layers-1) = Grid % zn(1:n_node_layers) - z_flat_bot
  dev_max = 0.0
  do n = 1, Grid % n_nodes
    lay = mod(n-1, n_node_layers)
    dev_max = max(dev_max, abs(Grid % zn(n) - z_flat_bot - s_lay(lay)))
  end do
  if(dev_max > MICRO * (z_flat_top - z_flat_bot)) then
    print '(a38,f12.4)', '# Largest deviation of levels [m]:    ', dev_max
    call Message % Error(83,                                             &
      'The node levels are not the same in all the columns of the '  //  &
      'grid.  The grid read by Convert must have a flat ground and ' //  &
      'the same vertical distribution of nodes everywhere.',             &
      file=__FILE__, line=__LINE__)
  end if

  ! First cell, it is defined in the .geo script now
  h_first = s_lay(1)
  cnt = count(col_rg(:) > h_first)
  print '(a38,f12.4)', '# Height of the first cell [m]:       ', h_first
  print '(a38,f12.4)', '# Safe height of the first cell [m]:  ', h_first_safe
  print '(a38,i9)',    '# Ground nodes rougher than that:     ', cnt
  if(h_first < h_first_safe) then
    print '(a)', ' # WARNING: The first cell is thinner than the safe height'
    print '(a)', ' # for this terrain.  With least squares gradients of'
    print '(a)', ' # pressure, simulation might diverge close to the roughest'
    print '(a)', ' # parts of the terrain.  Set GRADIENT_METHOD_FOR_PRESSURE'
    print '(a)', ' # to gauss_theorem in control file (stable on such grids),'
    print '(a)', ' # smooth the terrain (ROUGHNESS_ALLOWED in .ins file),'
    print '(a)', ' # thicken the first cell or refine the grid horizontally.'
  end if

  ! Thickness of the layer at the tallest roof
  dz_bld = h_first
  do lay = 1, n_node_layers - 1
    if(s_lay(lay) >= max_obstacle_height - MICRO) then
      dz_bld = s_lay(lay) - s_lay(lay-1)
      exit
    end if
  end do

  ! Zone in which roofs are flattened (and the layers of the flat grid
  ! kept).  It holds the tallest building and h_above over it.
  h_above = max(0.5 * max_obstacle_height, 5.0 * dz_bld)
  if(h_above_ins > -HUGE) h_above = h_above_ins
  k_zone = -1
  do lay = 1, n_node_layers - 2
    if(s_lay(lay) >= max_obstacle_height + h_above - MICRO) then
      k_zone = lay
      exit
    end if
  end do
  if(k_zone .eq. -1) then
    call Message % Error(85,                                              &
      'The zone of buildings (the tallest building and the height '   //  &
      'above it) does not fit below the top of the grid.  Raise the ' //  &
      'top in the .geo script, or lower HEIGHT_ABOVE_ROOFS.',             &
      file=__FILE__, line=__LINE__)
  end if
  h_zone = s_lay(k_zone)

  print '(a38,f12.4)', '# Layer thickness at tallest roof [m]:', dz_bld
  print '(a38,f12.4)', '# Height above the tallest roof [m]:  ', h_above
  print '(a38,f12.4)', '# Height of the zone of buildings [m]:', h_zone
  print '(a38,i9)',    '# Node layers in that zone:           ', k_zone
  print '(a38,i9)',    '# Node layers above it:               ',  &
                       n_node_layers - 1 - k_zone

  !-----------------------------------------------------------------!
  !                                                                 !
  !   Phase III: Find cells in buildings, chimneys and porosities   !
  !                                                                 !
  !-----------------------------------------------------------------!

  !--------------------------------!
  !   Treat them all in one loop   !
  !--------------------------------!
  n_buildings  = 0
  n_porosities = 0
  allocate(cell_in_building(Grid % n_cells));  cell_in_building(:) = 0
  allocate(cell_in_chimney (Grid % n_cells));  cell_in_chimney (:) = 0
  Assert(allocated(Grid % por));               Grid % por      (:) = 0
  allocate(bld_height(999));                   bld_height      (:) = 0.0

  do run = FOR_BUILDINGS, FOR_POROSITIES
  ! Browse through cells "on the ground"
  do c = 1, Grid % n_cells, n_cell_layers
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

            ! Nodes will be placed at heights s_lay above the ground
            ! (Phase IV).  Roof of a building goes to the level nearest
            ! to its nominal height (it is placed at the exact height in
            ! Phase V) and each building is at least one cell high.  For
            ! the rest, take the cells which are below the nominal height.
            n_lay = cu - c  ! number of cells below this one
            if(n_lay >= k_zone) goto 1
            if(run .eq. FOR_BUILDINGS) then
              bld_height(region) = height
              if(n_lay > 0 .and.  &
                 0.5 * (s_lay(n_lay) + s_lay(n_lay+1)) > height) goto 1
            else
              if(s_lay(n_lay+1) > height + MICRO) goto 1
            end if

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
  !   Phase IV: Map the grid onto the STL terrain                 !
  !                                                               !
  !   In each column, the nodes are placed at heights s_lay       !
  !   above the ground up to the top of the zone of buildings.    !
  !   Above, the layers of the flat grid are scaled uniformly to  !
  !   reach the flat top of the domain.  The scaling factor is    !
  !   found in each column separately, since columns have         !
  !   different heights.  On a flat ground the grid is unchanged. !
  !                                                               !
  !---------------------------------------------------------------!

  print *, '#=================================================================='
  print *, '# Terrain mapping: prescribed layers in the zone of buildings     '
  print *, '#------------------------------------------------------------------'
  print '(a38,2f12.4)', '# Flat-grid bottom and top:          ',  &
                        z_flat_bot, z_flat_top

  r_min =  HUGE
  r_max = -HUGE
  s_top = s_lay(n_node_layers-1)

  do col = 1, n_cols
    n        = (col-1) * n_node_layers + 1  ! ground node
    z_ground = col_zg(col)

    ! Space left above the zone of buildings
    l_up = z_flat_top - z_ground - h_zone
    if(l_up <= 0.0) then
      print '(a,3f14.4)', ' # Invalid terrain column (x,y,z): ', &
                       Grid % xn(n), Grid % yn(n), z_ground
      call Message % Error(86,                                         &
        'Terrain is too high for the domain top and the zone of '  //  &
        'buildings.  Raise the top boundary of the grid.',             &
        file=__FILE__, line=__LINE__)
    end if
    scale_up = l_up / (s_top - h_zone)
    r_min = min(r_min, scale_up)
    r_max = max(r_max, scale_up)

    ! Place the nodes
    do lay = 0, k_zone
      Grid % zn(n+lay) = z_ground + s_lay(lay)
    end do
    do lay = k_zone + 1, n_node_layers - 2
      Grid % zn(n+lay) = z_ground + h_zone + scale_up * (s_lay(lay) - h_zone)
    end do
    Grid % zn(n+n_node_layers-1) = z_flat_top  ! remove round-off
  end do

  print '(a38,2f12.4)', '# Layers above zone scaled, min/max: ', r_min, r_max
  if(r_min < 1.0 / MAX_SCALE .or. r_max > MAX_SCALE) then
    print '(a)', ' # WARNING: Layers above the zone of buildings are scaled'
    print '(a)', ' # a lot.  Consider moving the top of the domain, or adding'
    print '(a)', ' # layers in the vertical direction.'
  end if

  !---------------------------------------------------------------!
  !                                                               !
  !   Phase V: Flatten the roofs of the buildings                 !
  !                                                               !
  !   Each roof is pinned to a constant elevation, and the nodes  !
  !   in the columns of the building are redistributed between    !
  !   the ground, the roof(s) and the undisturbed grid above.     !
  !                                                               !
  !   The roof is placed at the nominal height of the building,   !
  !   above the average elevation of the terrain beneath, unless  !
  !   that leaves too little of the wall on the uphill side of    !
  !   the building (or sinks the roof in the terrain altogether,  !
  !   which leads to cells with negative volumes).  If that       !
  !   happens, the roof is raised.                                !
  !                                                               !
  !   Some building IDs can be absent (e.g. 001, 002, 035),       !
  !   so never assume every ID has its own ground cells.          !
  !---------------------------------------------------------------!

  n_cols = Grid % n_nodes / n_node_layers
  allocate(column_seen(n_cols));            column_seen(:)   = 0
  allocate(col_n_pins (n_cols));            col_n_pins (:)   = 0
  allocate(col_pin_lay(MAX_PINS, n_cols));  col_pin_lay(:,:) = 0
  allocate(col_pin_z  (MAX_PINS, n_cols));  col_pin_z  (:,:) = 0.0
  allocate(z_old(0:n_node_layers-1))

  print '(a)', ' # Flattening the roofs of the buildings'
  print '(a)', ' # bld. layer   nodes   ground min/avg/max [m]'  //  &
               '      nominal    meshed  roof ASL    raised'

  !------------------------------------------------------!
  !   Find roof elevation for each building and pin it   !
  !   in all the columns (z lines) the building holds    !
  !------------------------------------------------------!
  do b = 1, n_buildings

    ! Number of solid cells above the ground cell = roof node level
    top = -1
    consistent_roof = .true.
    do c = 1, Grid % n_cells, n_cell_layers
      if(cell_in_building(c) .ne. b) cycle
      col_top = 0
      do cu = c, c + n_cell_layers - 1
        if(cell_in_building(cu) .ne. b) exit
        col_top = col_top + 1
      end do
      if(top .eq. -1) top = col_top
      if(col_top .ne. top) consistent_roof = .false.
    end do

    ! A gap in building numbering does not denote an actual building
    if(top .eq. -1) cycle

    if(.not. consistent_roof .or. top < 1 .or. top >= n_node_layers-1) then
      print '(a,i5,a)', ' # Building ', b,  &
        ': roof levels inconsistent; leaving terrain-mapped roof unchanged'
      cycle
    end if

    ! Browse through all the columns of the building (also those at its
    ! perimeter, shared with the fluid or with another building) to find
    ! the extents of the ground beneath and the average roof elevation
    n_roof = 0
    zg_min = HUGE;  zg_max = -HUGE;  zg_avg = 0.0;  z_avg = 0.0
    do c = 1, Grid % n_cells, n_cell_layers
      if(cell_in_building(c) .ne. b) cycle
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        col = (Grid % cells_n(i_nod, c) - 1) / n_node_layers + 1
        if(column_seen(col) .ne. b) then
          column_seen(col) = b
          n = (col-1) * n_node_layers + 1  ! ground node
          n_roof = n_roof + 1
          zg_min = min(zg_min, Grid % zn(n))
          zg_max = max(zg_max, Grid % zn(n))
          zg_avg = zg_avg + Grid % zn(n)
          z_avg  = z_avg  + Grid % zn(n+top)
        end if
      end do
    end do
    zg_avg = zg_avg / real(n_roof)
    z_avg  = z_avg  / real(n_roof)
    h_meshed = z_avg - zg_avg  ! height of the building as it was meshed
    height   = bld_height(b)   ! nominal height of the building
    z_avg    = zg_avg + height

    ! Roof must stay well above the highest ground under the building
    z_roof = max(z_avg, zg_max + MIN_WALL_FRACTION * height)

    print '(a,i4,i6,i8,3f9.2,4f10.2)', ' # ', b, top, n_roof,  &
          zg_min, zg_avg, zg_max, height, h_meshed,  &
          z_roof, z_roof - z_avg

    ! Pin the roof in all the columns of the building
    do c = 1, Grid % n_cells, n_cell_layers
      if(cell_in_building(c) .ne. b) cycle
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        col = (Grid % cells_n(i_nod, c) - 1) / n_node_layers + 1
        if(column_seen(col) .ne. -b) then
          column_seen(col) = -b
          col_n_pins(col) = col_n_pins(col) + 1
          Assert(col_n_pins(col) .le. MAX_PINS)
          col_pin_lay(col_n_pins(col), col) = top
          col_pin_z  (col_n_pins(col), col) = z_roof
        end if
      end do
    end do

  end do  ! through buildings

  !------------------------------------------------------------!
  !   Redistribute nodes in the columns holding pinned roofs   !
  !------------------------------------------------------------!
  n_fixed   = 0
  ratio_min =  HUGE
  ratio_max = -HUGE
  do col = 1, n_cols
    if(col_n_pins(col) .eq. 0) cycle
    n = (col-1) * n_node_layers + 1  ! ground node
    z_old(0:n_node_layers-1) = Grid % zn(n:n+n_node_layers-1)

    ! Sort the pins by layers (there are just a few of them)
    do i = 2, col_n_pins(col)
      do j = i, 2, -1
        if(col_pin_lay(j, col) < col_pin_lay(j-1, col)) then
          lay = col_pin_lay(j, col)
          z_b = col_pin_z  (j, col)
          col_pin_lay(j, col) = col_pin_lay(j-1, col)
          col_pin_z  (j, col) = col_pin_z  (j-1, col)
          col_pin_lay(j-1, col) = lay
          col_pin_z  (j-1, col) = z_b
        end if
      end do
    end do

    ! Columns shared by buildings whose roofs are at the same layer; the
    ! roofs can't be both flat there, take the mean of the two elevations
    i = 1
    do while(i < col_n_pins(col))
      if(col_pin_lay(i, col) .eq. col_pin_lay(i+1, col)) then
        col_pin_z(i, col) = 0.5 * (col_pin_z(i, col) + col_pin_z(i+1, col))
        do j = i+1, col_n_pins(col)-1
          col_pin_lay(j, col) = col_pin_lay(j+1, col)
          col_pin_z  (j, col) = col_pin_z  (j+1, col)
        end do
        col_n_pins(col) = col_n_pins(col) - 1
      else
        i = i + 1
      end if
    end do

    ! Stretch (or squeeze) the layers between
    ! the ground and the roofs to fit the pins
    lay_a = 0
    z_a   = z_old(0)

    ! Keep the first cell at its prescribed height (it also neighbors
    ! the fluid at the perimeter of building) if the roof allows it
    lay_b = col_pin_lay(1, col)
    z_b   = col_pin_z  (1, col)
    if(lay_b >= 2) then
      if(z_b - z_old(1) >= MIN_DZ_RATIO * (z_old(lay_b) - z_old(1))) then
        lay_a = 1
        z_a   = z_old(1)
      end if
    end if

    do i = 1, col_n_pins(col)
      lay_b = col_pin_lay(i, col)
      z_b   = col_pin_z  (i, col)

      ! Roofs of adjoining buildings must not cross each other
      dz_min = MIN_DZ_RATIO * (z_old(lay_b) - z_old(lay_a))
      if(z_b < z_a + dz_min) then
        z_b = z_a + dz_min
        n_fixed = n_fixed + 1
      end if

      do lay = lay_a + 1, lay_b
        Grid % zn(n+lay) = z_a + (z_b - z_a)                     &
                               * (z_old(lay)   - z_old(lay_a))   &
                               / (z_old(lay_b) - z_old(lay_a))
      end do
      lay_a = lay_b
      z_a   = z_b
    end do

    ! Blend into the undisturbed grid above the highest roof.  The
    ! displacement of the roof decays linearly to the top of the zone
    ! of buildings, so all the layers in between change equally.
    lay_b = k_zone
    ratio = (z_old(lay_b) - z_a) / (z_old(lay_b) - z_old(lay_a))
    if(ratio <= 0.0) then
      print '(a,2f14.4)', ' # Roof above the zone of buildings at (x,y): ', &
                          Grid % xn(n), Grid % yn(n)
      call Message % Error(89,                                           &
        'A roof ended above the zone of buildings.  Increase the ' //    &
        'height of the zone above the tallest building.',                &
        file=__FILE__, line=__LINE__)
    end if
    ratio_min = min(ratio_min, ratio)
    ratio_max = max(ratio_max, ratio)
    do lay = lay_a + 1, lay_b - 1
      Grid % zn(n+lay) = z_old(lay_b) - (z_old(lay_b) - z_old(lay)) * ratio
    end do

  end do  ! through columns

  print '(a38,2f12.4)', '# Layers above roofs scaled, min/max:',  &
                        ratio_min, ratio_max
  if(ratio_min < MIN_DZ_RATIO .or. ratio_max > 1.0 / MIN_DZ_RATIO) then
    print '(a)', ' # WARNING: Layers above the roofs are deformed a lot.'
    print '(a)', ' # Consider increasing the zone above the tallest building.'
  end if

  if(n_fixed > 0) then
    print '(a,i6,a)', ' # WARNING: roofs of adjoining buildings were ' //  &
                      'too close in ', n_fixed, ' columns; moved apart'
  end if

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

  !------------------------------------------------------------------------!
  !                                                                        !
  !   Phase VII: Insert new boundary faces on building and chimney walls   !
  !                                                                        !
  !------------------------------------------------------------------------!
  if(n_buildings .gt. 0) then

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)
      if(c2 .gt. 0) then

        c1_on = cell_in_building(c1) .eq. 0 ! .and. cell_in_chimney(c1) .eq. 0
        c2_on = cell_in_building(c2) .eq. 0 ! .and. cell_in_chimney(c2) .eq. 0

        if(      c1_on .and.       c2_on) cycle
        if(.not. c1_on .and. .not. c2_on) cycle

        if(c1_on .and. .not. c2_on) c = c1
        if(c2_on .and. .not. c1_on) c = c2

        s_nod(1:4) = Grid % faces_n(1:4, s)
        call Sort % Int_Array(s_nod)

        if(Grid % cells_n_nodes(c) .eq. 4) fn = TET
        if(Grid % cells_n_nodes(c) .eq. 5) fn = PYR
        if(Grid % cells_n_nodes(c) .eq. 6) fn = WED
        if(Grid % cells_n_nodes(c) .eq. 8) fn = HEX

        do dir = 1, 6
          if(Grid % cells_bnd_region(dir, c) .eq. 0) then

            n_f_nod    = 0
            f_nod(1:4) = 0

            ! Fill up face nodes
            do i_nod = 1, 4
              if(fn(dir, i_nod) > 0) then
                f_nod(i_nod) = Grid % cells_n(fn(dir, i_nod), c)
                n_f_nod   = n_f_nod + 1
              end if
            end do
            call Sort % Int_Array(f_nod)

            ! Check all nodes from both faves
            if (all(s_nod == f_nod)) then
              Grid % cells_bnd_region(dir, c) = 1
              exit
            end if

          end if  ! cells_bnd_region(dir, c) .eq. 0
        end do    ! dir
      end if      ! c2 .gt .0
    end do        ! n_buildings .gt. 0

  end if

  do bc = Boundary_Regions()
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name(1:7) .eq. 'CHIMNEY') then

      read(bc_name(13:15), *) region  ! this is, in essence, building number
                                      ! (4 years l8r: is it chimney number?)
      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1, s)
        c2 = Grid % faces_c(2, s)
        if(c2 .gt. 0) then

          c1_on = cell_in_chimney(c1) .eq. 0
          c2_on = cell_in_chimney(c2) .eq. 0

          if(      c1_on .and.       c2_on) cycle
          if(.not. c1_on .and. .not. c2_on) cycle

          if(c1_on .and. .not. c2_on) c = c1
          if(c2_on .and. .not. c1_on) c = c2

          s_nod(1:4) = Grid % faces_n(1:4, s)
          call Sort % Int_Array(s_nod)

          if(Grid % cells_n_nodes(c) .eq. 4) fn = TET
          if(Grid % cells_n_nodes(c) .eq. 5) fn = PYR
          if(Grid % cells_n_nodes(c) .eq. 6) fn = WED
          if(Grid % cells_n_nodes(c) .eq. 8) fn = HEX

          do dir = 1, 6
            if(Grid % cells_bnd_region(dir, c) .eq. 0) then

              n_f_nod    = 0
              f_nod(1:4) = 0

              ! Fill up face nodes
              do i_nod = 1, 4
                if(fn(dir, i_nod) > 0) then
                  f_nod(i_nod) = Grid % cells_n(fn(dir, i_nod), c)
                  n_f_nod   = n_f_nod + 1
                end if
              end do
              call Sort % Int_Array(f_nod)

              ! Check all nodes from both faves
              if (all(s_nod == f_nod)) then
                Grid % cells_bnd_region(dir, c) = bc
                exit
              end if

            end if  ! cells_bnd_region(dir, c) .eq. 0
          end do    ! dir
        end if      ! c2 .gt. 0
      end do        ! s

    end if  ! chimney

  end do  ! boundary regions

  !-------------------------------------------------------------!
  !                                                             !
  !   Phase VIII: Store only cells which are not in buildings   !
  !                                                             !
  !-------------------------------------------------------------!

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

  !----------------------------------------------------!
  !                                                    !
  !   Phase IX: Smooth the grid lines in z direction   !
  !                                                    !
  !----------------------------------------------------!

  print '(a)', ' # Smoothing the grid out a little bit'

  ! This is a wee-bit memory intensive, better be on the safe side
  deallocate(cell_in_building)
  deallocate(cell_in_chimney)
  deallocate(node_on_chimney)
  deallocate(column_seen)
  deallocate(col_pin_lay)
  deallocate(col_pin_z)
  deallocate(z_old)
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

  ! Keep the first cell at the prescribed height above the ground
  do n = 1, Grid % n_nodes
    if(mod(n-1, n_node_layers) .eq. 1) nodes_n_nodes(n) = 0
  end do

  deallocate(col_n_pins)

  !---------------------------!
  !   Find max z coordinate   !
  !---------------------------!
  z_max = -HUGE
  do n = 1, Grid % n_nodes
    z_max = max(z_max, Grid % zn(n))
  end do

  !-------------------------------------------------------------------!
  !   Smooth in the z direction a bit.                                !
  !                                                                   !
  !   This damps the roughness of the terrain with the distance from  !
  !   the ground, and the steps introduced by flattening of the roofs !
  !   with the distance from the buildings.  On a rough terrain, the  !
  !   smoothing pulls the nodes close to the ground towards their     !
  !   neighbors, which squeezes (and, if the layers are thin, also    !
  !   inverts) the cells.  Hence, no layer is allowed to get thinner  !
  !   than a fraction of the thickness it had before the smoothing.   !
  !-------------------------------------------------------------------!
  ! Number of iterations: Z_SMOOTH_ITERS in .ins file, environment variable
  ! CITY_Z_SMOOTH_ITERS, or ITERS (0 keeps the z coordinates of Phase V)
  print '(a,i4)', ' # Iterations of smoothing in z direction: ', n_smooth
  allocate(z_ref(Grid % n_nodes))
  z_ref(1:Grid % n_nodes) = Grid % zn(1:Grid % n_nodes)

  do iter = 1, n_smooth

    ! Compute averaged z coordinate for each node
    do n = 1, Grid % n_nodes
      if(nodes_n_nodes(n) > 0 .and. Grid % zn(n) < z_max) then

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
      if(nodes_n_nodes(n) > 0 .and. Grid % zn(n) < z_max) then

        ! No matter how small, some under-relaxation is needed here
        z_avg = 0.1 * Grid % zn(n) + 0.9 * z_new(n)

        ! Nodes n-1 and n+1 are below and above node n, in the same column
        ! (Nodes on the ground and at the top are not smoothed, see above.)
        z_a = Grid % zn(n-1) + MIN_SMOOTH_DZ * (z_ref(n)   - z_ref(n-1))
        z_b = Grid % zn(n+1) - MIN_SMOOTH_DZ * (z_ref(n+1) - z_ref(n))
        Grid % zn(n) = min(max(z_avg, z_a), z_b)
      end if
    end do  ! through nodes
  end do    ! iterations

  call Profiler % Stop('Insert_Buildings')

  end subroutine
