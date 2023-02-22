!==============================================================================!
  subroutine Insert_Buildings(Convert, Grid)
!------------------------------------------------------------------------------!
!   This soubroutine was developed to import mesh with buildings, and place    !
!   them both on a ground.  This subroutine is not a part of the standard      !
!   T-Flows distribution, it was developed specifically for the project on     !
!   Smart Cities.  It is rather long and it will probably stay like that.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert
  type(Grid_Type)     :: Grid
!-----------------------------------[Locals]-----------------------------------!

  ! Facet type
  type Facet_Type
    real :: x(3), y(3), z(3)  ! vertex's coordinates
    real :: area_z
  end type

  ! Ground type
  type Ground_Type
    integer                       :: n_facets
    type(Facet_Type), allocatable :: facet(:)
  end type

  ! Variables for phase I
  integer           :: cnt, fu, v
  real              :: x_min, x_max, y_min, y_max, z_min, z_max, area
  character(len=SL) :: name_in
  type(Ground_Type) :: ground

  ! Variables for phase II (in addition to exisiting ones)
  integer              :: n, f
  real                 :: area_1, area_2, area_3, z, z_scale, z_from_top
  integer, allocatable :: ground_cell(:)
  character(len=SL)    :: bc_name

  ! Variables for phase III (in addition to exisiting ones)
  logical :: cells_in_columns
  integer :: run, c, cg, dir, n_ground_cells, near_cg
  real    :: dis2, dis2_min

  ! Variables for phase IV (in addition to exisiting ones)
  integer, allocatable :: i_work_1(:)
  integer, allocatable :: i_work_2(:,:)
  integer, allocatable :: i_work_3(:,:)
  integer, allocatable :: new_c(:)
  integer, allocatable :: old_c(:)
  integer, allocatable :: criteria(:,:)

  ! Variables for phase V (in addition to exisiting ones)
  integer              :: cu                   ! cell upper
  integer              :: n_hor_layers, ni
  logical, allocatable :: cell_in_building(:)
  logical, allocatable :: node_on_building(:)
  integer, allocatable :: cell_in_chimney(:)
  integer, allocatable :: node_on_chimney(:)
  logical              :: buildings_exist
  integer              :: n_chimneys, chimney
  real                 :: height               ! building height

  ! Variables for phase VI (in addition to exisiting ones)
  integer :: bc, ground_bc

  ! No new variables for phase VII
  ! Variables for phase VIII (in addition to exisiting ones)
  integer :: fn(6,4), n_f_nod, f_nod(4)
!==============================================================================!

  !----------------------------------------------------------!
  !                                                          !
  !   Phase I: Read the ground definition from an STL file   !
  !                                                          !
  !----------------------------------------------------------!

  print *, '#==============================================================='
  print *, '# Enter the name of the ground file (with ext.):'
  print *, '#---------------------------------------------------------------'
  read(*,*) name_in

  call File % Open_For_Reading_Ascii(name_in, fu)

  !--------------------------!
  !   Count all the facets   !
  !--------------------------!
  rewind(fu)
  cnt = 0
  do
    call File % Read_Line(fu)
    if(Line % tokens(1) .eq. 'endsolid') exit
    if(Line % tokens(1) .eq. 'facet') cnt = cnt + 1
  end do
  print '(a38,i9)', '# Number of facets on the ground:    ', cnt

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  ground % n_facets = cnt
  allocate(ground % facet(ground % n_facets))

  !---------------------------------!
  !   Read all vertex coordinates   !
  !---------------------------------!
  x_max = -HUGE;  y_max = -HUGE;  z_max = -HUGE
  x_min = +HUGE;  y_min = +HUGE;  z_min = +HUGE
  area = 0.0
  rewind(fu)
  cnt = 0
  do
    call File % Read_Line(fu)
    if(Line % tokens(1) .eq. 'endsolid') exit
    if(Line % tokens(1) .eq. 'facet') then
      cnt = cnt + 1
      call File % Read_Line(fu)                ! 'outer loop'
      do v = 1, 3
        call File % Read_Line(fu)              ! 'vertex 1, 2 and 3'
        read(Line % tokens(2), *) ground % facet(cnt) % x(v)
        read(Line % tokens(3), *) ground % facet(cnt) % y(v)
        read(Line % tokens(4), *) ground % facet(cnt) % z(v)
        x_max = max(x_max, ground % facet(cnt) % x(v))
        x_min = min(x_min, ground % facet(cnt) % x(v))
        y_max = max(y_max, ground % facet(cnt) % y(v))
        y_min = min(y_min, ground % facet(cnt) % y(v))
        z_max = max(z_max, ground % facet(cnt) % z(v))
        z_min = min(z_min, ground % facet(cnt) % z(v))
      end do
      call File % Read_Line(fu)                ! 'endloop'
      ground % facet(cnt) % area_z =                     &
        Convert % Triangle_Area_Z(ground % facet(cnt) % x(1),  &
                                  ground % facet(cnt) % y(1),  &
                                  ground % facet(cnt) % x(2),  &
                                  ground % facet(cnt) % y(2),  &
                                  ground % facet(cnt) % x(3),  &
                                  ground % facet(cnt) % y(3))
      area = area + ground % facet(cnt) % area_z
    end if
  end do
  print *, '# Read all ground facets!'
  print '(a38,6f12.4)', '# Bounding box: (corner 1)           ',  &
                                                       x_min, y_min, z_min
  print '(a38,6f12.4)', '#               (corner 2)           ',  &
                                                       x_max, y_max, z_max
  print '(a38, f12.4)', '# Ground area from bounding box:     ',  &
                                          (x_max - x_min) * (y_max - y_min)
  print '(a38, f12.4)', '# Ground area as a sum of facets:    ', area

  !---------------------------------------------------------------!
  !                                                               !
  !   Phase II: Place grid on the terrain defind in an STL file   !
  !                                                               !
  !---------------------------------------------------------------!

  print *, '#=================================================================='
  print *, '# Placing grid on the ground.  This may take a few minutes         '
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
  do n = 1, Grid % n_nodes

    ! Print progress on the screen
    if(mod(n, (Grid % n_nodes / 20) ) .eq. 0) then
      print '(a2, f5.0, a14)',   &
            ' #', (100. * n / (1.0*(Grid % n_nodes))), ' % complete...'
    end if ! each 5%

    do f = 1, ground % n_facets

      area = ground % facet(f) % area_z  ! area has other meaning in this phase

      area_1 = Convert % Triangle_Area_Z(Grid % xn(n), Grid % yn(n),    &
         ground % facet(f) % x(2), ground % facet(f) % y(2),            &
         ground % facet(f) % x(3), ground % facet(f) % y(3))

      area_2 = Convert % Triangle_Area_Z(Grid % xn(n), Grid % yn(n),    &
         ground % facet(f) % x(3), ground % facet(f) % y(3),            &
         ground % facet(f) % x(1), ground % facet(f) % y(1))

      area_3 = Convert % Triangle_Area_Z(Grid % xn(n), Grid % yn(n),    &
         ground % facet(f) % x(1), ground % facet(f) % y(1),            &
         ground % facet(f) % x(2), ground % facet(f) % y(2))

      if((area_1+area_2+area_3 - NANO) < area) then
        z = (area_1 / area) * ground % facet(f) % z(1)  &
          + (area_2 / area) * ground % facet(f) % z(2)  &
          + (area_3 / area) * ground % facet(f) % z(3)
        z_from_top   = z_max - Grid % zn(n)
        z_scale      = (z_max - z) / (z_max - z_min)
        Grid % zn(n) = z_max - z_from_top * z_scale
      end if
    end do
  end do

  !--------------------------------------------------!
  !                                                  !
  !   Phase III: Align horizontal cell coordinates   !
  !                                                  !
  !--------------------------------------------------!

  !---------------------------------------------!
  !   Calculate cell centers because you will   !
  !     be sorting them by height later on      !
  !---------------------------------------------!
  call Grid % Calculate_Cell_Centers()

  !--------------------------------------------------!
  !   Check if cells are already stored in columns   !
  !--------------------------------------------------!
  cells_in_columns = .true.
  do c = 1, Grid % n_cells
    if(Grid % zc(c+1) < Grid % zc(c)) then
      n_hor_layers = c
      exit
    end if
  end do
  print '(a38,i9)', '# Tentative number of layers:        ', n_hor_layers
  if(mod(Grid % n_cells, n_hor_layers) .ne. 0) cells_in_columns = .false.
  do c = 1, Grid % n_cells, n_hor_layers
    if(Grid % zc(c+1) < Grid % zc(c))          cells_in_columns = .false.
  end do
  do c = n_hor_layers, Grid % n_cells, n_hor_layers
    if(Grid % zc(c-1) > Grid % zc(c))          cells_in_columns = .false.
  end do
  if(cells_in_columns) goto 1

  !----------------------------------------!
  !   Find and store cells on the ground   !
  !----------------------------------------!
  do run = 1, 2
    n_ground_cells = 0
    do c = 1, Grid % n_cells
      do dir = 1, 6
        if(Grid % cells_bnd_region(dir, c) .ne. 0) then
          bc_name = trim(Grid % region % name(Grid % cells_bnd_region(dir, c)))
          call String % To_Upper_Case(bc_name)
          if(bc_name(1:8) .eq. 'BUILDING' .or.  &
             bc_name      .eq. 'GROUND') then

            n_ground_cells = n_ground_cells + 1
            if(run .eq. 2) ground_cell(n_ground_cells) = c
          end if
        end if
      end do
    end do
    if(run .eq. 1) then
      allocate(ground_cell(n_ground_cells)); ground_cell(:) = 0
    end if
  end do
  print '(a38,i9)', '# Number of cells on the ground:     ', n_ground_cells
  print '(a38,i9)', '# Number of horizontal layers (1)    ', Grid % n_cells  &
                                                           / n_ground_cells

  !---------------------------------------------------------------------!
  !   For each cell, find the nearest on the ground and align with it   !
  !---------------------------------------------------------------------!

  print *, '#=================================================================='
  print *, '# Aligning cell coordinates.  This may take a few minutes          '
  print *, '#------------------------------------------------------------------'
  do c = 1, Grid % n_cells
    dis2_min = HUGE

    ! Print progress on the screen
    if(mod(c, (Grid % n_cells / 20) ) .eq. 0) then
      print '(a2, f5.0, a14)',   &
            ' #', (100. * c / (1.0*(Grid % n_cells))), ' % complete...'
    end if ! each 5%

    do n = 1, n_ground_cells
      cg = ground_cell(n)     ! real ground cell number

      ! Compute planar distance
      dis2 = (  (Grid % xc(c) - Grid % xc(cg)) ** 2  &
              + (Grid % yc(c) - Grid % yc(cg)) ** 2)

      ! Store ground cell if nearest
      if(dis2 < dis2_min) then
        dis2_min = dis2
        near_cg = cg
      end if
    end do

    ! Simply take planar coordinates from the nearest wall cell
    Grid % xc(c) = Grid % xc(near_cg)
    Grid % yc(c) = Grid % yc(near_cg)
  end do

  !----------------------------------------------------!
  !                                                    !
  !   Phase IV: Sort cells in columns in z direction   !
  !                                                    !
  !----------------------------------------------------!

  allocate(i_work_1(   Grid % n_cells));     i_work_1(:)   = 0
  allocate(i_work_2(8, Grid % n_cells));     i_work_2(:,:) = 0
  allocate(i_work_3(6, Grid % n_cells));     i_work_3(:,:) = 0
  allocate(new_c   (   Grid % n_cells));     new_c(:)      = 0
  allocate(old_c   (   Grid % n_cells));     old_c(:)      = 0
  allocate(criteria(   Grid % n_cells, 3));  criteria(:,:) = 0

  !------------------------!
  !   Store cells' nodes   !
  !------------------------!
  do c = 1, Grid % n_cells
    i_work_1(c)      = Grid % cells_n_nodes(c)
    i_work_2(1:8, c) = Grid % cells_n(1:8, c)
    i_work_3(1:6, c) = Grid % cells_bnd_region(1:6, c)
  end do

  !--------------------------!
  !   Set sorting criteria   !
  !--------------------------!
  do c = 1, Grid % n_cells
    criteria(c, 1) = nint(Grid % xc(c) * MEGA)
    criteria(c, 2) = nint(Grid % yc(c) * MEGA)
    criteria(c, 3) = nint(Grid % zc(c) * MEGA)
    old_c(c)       = c
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort % Three_Int_Carry_Int(criteria(1:Grid % n_cells, 1),  &
                                  criteria(1:Grid % n_cells, 2),  &
                                  criteria(1:Grid % n_cells, 3),  &
                                  old_c   (1:Grid % n_cells))

  ! This was a bit of a bluff but it worked
  do c = 1, Grid % n_cells
    new_c(old_c(c)) = c
  end do

  !-----------------------------------------------!
  !   Do the sorting of data pertinent to cells   !
  !-----------------------------------------------!
  do c = 1, Grid % n_cells
    Grid % cells_n_nodes        (new_c(c)) = i_work_1(c)
    Grid % cells_n         (1:8, new_c(c)) = i_work_2(1:8, c)
    Grid % cells_bnd_region(1:6, new_c(c)) = i_work_3(1:6, c)
  end do
  call Sort % Real_By_Index(Grid % n_cells, Grid % xc(1), new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % yc(1), new_c(1))
  call Sort % Real_By_Index(Grid % n_cells, Grid % zc(1), new_c(1))

  deallocate(i_work_1)
  deallocate(i_work_2)
  deallocate(i_work_3)
  deallocate(new_c   )
  deallocate(old_c   )
  deallocate(criteria)

1 continue

  !--------------------------------------!
  !                                      !
  !   Phase V: Find cells in buildings   !
  !                                      !
  !--------------------------------------!

  ! Estimate number of horizontal layers again (to check if sorting worked)
  do c = 1, Grid % n_cells
    if(Grid % zc(c+1) < Grid % zc(c)) then
      n_hor_layers = c
      exit
    end if
  end do
  print '(a38,i9)', '# Number of horizontal layers (2):   ', n_hor_layers

  !-----------------------------!
  !   Find cells in buildings   !
  !-----------------------------!
  buildings_exist = .false.
  allocate(cell_in_building(Grid % n_cells))
  cell_in_building(:) = .false.
  do c = 1, Grid % n_cells
    do dir = 1, 6
      if(Grid % cells_bnd_region(dir, c) .ne. 0) then
        bc_name = trim(Grid % region % name(Grid % cells_bnd_region(dir, c)))
        call String % To_Upper_Case(bc_name)
        if(bc_name(1:8) .eq. 'BUILDING') then
          read(bc_name(10:12), *) height
!          height = height / 1000.0
          do cu = c, c + n_hor_layers - 1
            if(Grid % zc(cu) - height >= Grid % zc(c)) exit
            cell_in_building(cu) = .true.        ! mark cell as in building
            buildings_exist = .true.             ! at least one building found
          end do
        end if
      end if
    end do
  end do

  !----------------------------!
  !   Find cells in chimneys   !
  !----------------------------!
  n_chimneys = 0
  allocate(cell_in_chimney(Grid % n_cells))
  cell_in_chimney(:) = 0
  do bc = 1, Grid % n_bnd_regions
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name(1:7) .eq. 'CHIMNEY') then

      read(bc_name(9:11), *)  height
      read(bc_name(13:15), *) chimney

      n_chimneys = n_chimneys + 1

      do c = 1, Grid % n_cells
        do dir = 1, 6
          if(Grid % cells_bnd_region(dir, c) .eq. bc) then
            do cu = c, c + n_hor_layers - 1
              if(Grid % zc(cu) - height >= Grid % zc(c)) exit
              cell_in_chimney(cu) = chimney  ! mark cell as in this chimney
            end do
          end if
        end do
      end do

    end if
  end do

  ! Count the number of cells in buildings (this is for checking only)
  cnt = 0
  do c = 1, Grid % n_cells
    if(cell_in_building(c)) cnt = cnt + 1
  end do
  print '(a38,i9)', '# Number of cells in buildings:      ', cnt

  !-----------------------------!
  !   Mark nodes on buildings   !
  !-----------------------------!
  allocate(node_on_building(Grid % n_nodes))
  node_on_building(:) = .false.
  do c = 1, Grid % n_cells
    if(cell_in_building(c)) then
      do ni = 1, Grid % cells_n_nodes(c)  ! mark node as on building
        n = Grid % cells_n(ni, c)
        node_on_building(n) = .true.
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
      do ni = 1, Grid % cells_n_nodes(c)  ! mark node as on chimney
        n = Grid % cells_n(ni, c)
        node_on_chimney(n) = cell_in_chimney(c)
      end do
    end if
  end do

  !-----------------------------------------------------------!
  !                                                           !
  !   Phase VI: Manage boundary condition names and numbers   !
  !                                                           !
  !-----------------------------------------------------------!

  ! Find ground b.c. number (only used to eliminate BUILDING_000)
  do bc = 1, Grid % n_bnd_regions
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name .eq. 'GROUND') then
      ground_bc = bc
    end if
  end do

  ! Eliminate building_000 b.c.
  do bc = 1, Grid % n_bnd_regions
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name .eq. 'BUILDING_000') then
      do c = 1, Grid % n_cells
        do dir = 1, 6
          if(Grid % cells_bnd_region(dir, c) .eq. bc) then
            Grid % cells_bnd_region(dir, c) = ground_bc
          end if
        end do
      end do
    end if
  end do

  cnt   = 0
  do bc = 1, Grid % n_bnd_regions
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name(1:8) .ne. 'BUILDING') then
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
  if(buildings_exist) then
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

  allocate(i_work_1(   Grid % n_cells));     i_work_1(:)   = 0
  allocate(i_work_2(8, Grid % n_cells));     i_work_2(:,:) = 0
  allocate(i_work_3(6, Grid % n_cells));     i_work_3(:,:) = 0
  allocate(new_c   (   Grid % n_cells));     new_c(:)      = 0
  allocate(old_c   (   Grid % n_cells));     old_c(:)      = 0
  allocate(criteria(   Grid % n_cells, 3));  criteria(:,:) = 0

  !-----------------------------------------------------!
  !   Store cells' nodes and boundary conditons again   !
  !-----------------------------------------------------!
  do c = 1, Grid % n_cells
    i_work_1(c)      = Grid % cells_n_nodes(c)
    i_work_2(1:8, c) = Grid % cells_n(1:8, c)
    i_work_3(1:6, c) = Grid % cells_bnd_region(1:6, c)
  end do

  !------------------------------------------------------!
  !   Renumber the cells without the cells in building   !
  !------------------------------------------------------!
  new_c(:) = 0
  cnt      = 0
  do c = 1, Grid % n_cells
    if(.not. cell_in_building(c)            &
       .and. cell_in_chimney(c) == 0) then
      cnt = cnt + 1
      new_c(c) = cnt
    end if
  end do
  print '(a38,i9)', '# Old number of cells:               ', Grid % n_cells
  print '(a38,i9)', '# Number of cells without buildings: ', cnt
  print '(a38,i9)', '# Number of excluded cells:          ', Grid % n_cells-cnt

  !--------------------------------------!
  !   Information on cell connectivity   !
  !--------------------------------------!
  do c = 1, Grid % n_cells
    if(new_c(c) .ne. 0) then
      Grid % cells_n_nodes       (new_c(c)) = i_work_1(c)
      Grid % cells_n        (1:8, new_c(c)) = i_work_2(1:8, c)
      Grid % cells_bnd_region(1:6, new_c(c)) = i_work_3(1:6, c)
    end if
  end do
  Grid % n_cells = cnt

  deallocate(i_work_1)
  deallocate(i_work_2)
  deallocate(i_work_3)
  deallocate(new_c   )
  deallocate(old_c   )
  deallocate(criteria)

  !-------------------------------------------------------------------------!
  !                                                                         !
  !   Phase VIII: Insert new boundary faces on building and chimney walls   !
  !                                                                         !
  !-------------------------------------------------------------------------!
  if(buildings_exist) then

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
          do ni = 1, 4
            if(fn(dir, ni) > 0) then
              f_nod(ni) = Grid % cells_n(fn(dir, ni), c)
              n_f_nod   = n_f_nod + 1
            end if
          end do

          if( n_f_nod > 0 ) then

            ! Qadrilateral face
            if(f_nod(4) > 0) then
              if( node_on_building(f_nod(1)) .and.  &
                  node_on_building(f_nod(2)) .and.  &
                  node_on_building(f_nod(3)) .and.  &
                  node_on_building(f_nod(4)) ) then
                Grid % cells_bnd_region(dir, c) = 1
              end if
            else
              if( node_on_building(f_nod(1)) .and.  &
                  node_on_building(f_nod(2)) .and.  &
                  node_on_building(f_nod(3)) ) then
                Grid % cells_bnd_region(dir, c) = 1
              end if
            end if
          end if
        end if
      end do
    end do

  end if

  do bc = 1, Grid % n_bnd_regions
    bc_name = trim(Grid % region % name(bc))
    call String % To_Upper_Case(bc_name)
    if(bc_name(1:7) .eq. 'CHIMNEY') then

      read(bc_name(13:15), *) chimney  ! this is, in essence, building number

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
            do ni = 1, 4
              if(fn(dir, ni) > 0) then
                f_nod(ni) = Grid % cells_n(fn(dir, ni), c)
                n_f_nod   = n_f_nod + 1
              end if
            end do

            if( n_f_nod > 0 ) then

              ! Qadrilateral face
              if(f_nod(4) > 0) then
                if( node_on_chimney(f_nod(1)) .eq. chimney .and.  &
                    node_on_chimney(f_nod(2)) .eq. chimney .and.  &
                    node_on_chimney(f_nod(3)) .eq. chimney .and.  &
                    node_on_chimney(f_nod(4)) .eq. chimney ) then
                  Grid % cells_bnd_region(dir, c) = bc
                end if
              else
                if( node_on_chimney(f_nod(1)) .eq. chimney .and.  &
                    node_on_chimney(f_nod(2)) .eq. chimney .and.  &
                    node_on_chimney(f_nod(3)) .eq. chimney ) then
                  Grid % cells_bnd_region(dir, c) = bc
                end if
              end if
            end if
          end if
        end do
      end do

    end if

  end do  ! chimney

  end subroutine
