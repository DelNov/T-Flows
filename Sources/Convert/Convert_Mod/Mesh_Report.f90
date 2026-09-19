!==============================================================================!
  subroutine Mesh_Report(Convert, Grid)
!------------------------------------------------------------------------------!
!>  Compiles a report on the converted mesh: size, boundary regions, wall
!>  distance, quality metrics (non-orthogonality, skewness, aspect ratio and
!>  volume ratio of neighbouring cells) and a regularity check of all cells
!>  (positive volumes, closed and outward oriented cells, consistency of cell
!>  shapes, unused nodes, disconnected domains).  The report is printed on
!>  the screen and stored in the file "<name>-mesh_report.txt".
!------------------------------------------------------------------------------!
!   Notes:                                                                     !
!                                                                              !
!   * It must be called after the geometry has been calculated, cells and      !
!     faces sorted and Grid % Find_Cells_Faces has been called.                !
!   * Wall distance is reported only if it was calculated.                     !
!   * Quality of periodic faces (those with shadows) is not assessed, because  !
!     the two cells they connect are not neighbours in the physical space.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
!---------------------------------[Parameters]---------------------------------!
  real, parameter :: NON_ORTH_WARN = 70.0     ! degrees
  real, parameter :: NON_ORTH_BAD  = 85.0
  real, parameter :: SKEW_WARN     = 0.8
  real, parameter :: SKEW_BAD      = 0.95
  real, parameter :: ASPECT_WARN   = 1000.0
  real, parameter :: VOL_RAT_WARN  = 20.0
  real, parameter :: CLOSURE_TOL   = 1.0e-6
  integer, parameter :: MAX_MSG    = 40
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s, sf, i_fac, n, b, i, fu, nb, nc
  integer :: n_int_faces, n_bnd_faces, n_per_faces, n_comp, n_seen, head, tail
  integer :: n_warn, n_err
  integer :: n_vol_bad, n_open, n_inward, n_few_faces, n_shape, n_bad_node
  integer :: n_zero_area, n_bad_link, n_bad_region, n_inverted, n_unused
  integer :: f_vol_bad, f_open, f_inward, f_few, f_shape, f_zero_area
  integer :: n_no, n_sk, n_ar, n_vr, w_no, w_sk, w_ar, w_vr, min_comp
  integer :: hist_no(6), hist_sk(6)
  real    :: xc1, yc1, zc1, xd, yd, zd, sx, sy, sz, ss, dd, ds, t, cs
  real    :: x_i, y_i, z_i, angle, skew, aspect, ratio, vol1, vol2, dn
  real    :: sum_x, sum_y, sum_z, sum_s, dist, d_min, d_max, s_sign
  real    :: no_min, no_max, no_sum, sk_min, sk_max, sk_sum
  real    :: ar_max, ar_sum, vr_max, vr_sum
  real    :: no_max_at(3), sk_max_at(3), ar_max_at(3), vr_max_at(3)
  real    :: wd_min, wd_min_at(3), min_vol, max_vol, tot_vol, min_vol_at(3)
  real    :: xmin, ymin, zmin, xmax, ymax, zmax, bin_edge_no(5), bin_edge_sk(5)
  integer :: n_wd, n_wd_unset
  logical :: bad_cell
  character(SL)             :: name_out
  character(256)            :: line, tmp
  character(256)            :: msg(MAX_MSG)
  logical,       allocatable :: node_used(:), visited(:)
  integer,       allocatable :: queue(:), b_count(:)
  real,          allocatable :: b_area(:), b_dmin(:), b_dmax(:), b_dsum(:)
  real,          allocatable :: cell_flag(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  call Profiler % Start('Mesh_Report')

  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  n_warn = 0
  n_err  = 0

  bin_edge_no = (/20.0, 40.0, 60.0, 70.0, 80.0/)
  bin_edge_sk = (/0.1, 0.25, 0.5, 0.8, 0.95/)

  allocate(cell_flag(-nb:nc));  cell_flag(:) = 0.0
  allocate(b_count(Grid % n_bnd_regions));  b_count(:) = 0
  allocate(b_area (Grid % n_bnd_regions));  b_area (:) = 0.0
  allocate(b_dmin (Grid % n_bnd_regions));  b_dmin (:) = +HUGE
  allocate(b_dmax (Grid % n_bnd_regions));  b_dmax (:) = -HUGE
  allocate(b_dsum (Grid % n_bnd_regions));  b_dsum (:) = 0.0

  !-------------------------------------!
  !   Open the file for storing report  !
  !-------------------------------------!
  call File % Set_Name(name_out, appendix='-mesh_report', extension='.txt')
  call File % Open_For_Writing_Ascii(name_out, fu)

  !-----------------------------------------!
  !                                         !
  !   Pass through faces: quality metrics   !
  !                                         !
  !-----------------------------------------!
  n_int_faces = 0;  n_bnd_faces = 0;  n_per_faces = 0
  n_zero_area = 0;  n_bad_link  = 0;  n_bad_region = 0;  n_inverted = 0
  f_zero_area = 0
  n_no = 0;  no_min = +HUGE;  no_max = -HUGE;  no_sum = 0.0;  w_no = 0
  n_sk = 0;  sk_min = +HUGE;  sk_max = -HUGE;  sk_sum = 0.0;  w_sk = 0
  n_vr = 0;  vr_max = 1.0;    vr_sum = 0.0;    w_vr = 0
  no_max_at = 0.0;  sk_max_at = 0.0;  vr_max_at = 0.0
  hist_no = 0;  hist_sk = 0

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c1 .le. 0 .or. c1 .gt. nc .or. c2 .eq. 0 .or. c2 .gt. nc) then
      n_bad_link = n_bad_link + 1
      cycle
    end if
    if(c1 .eq. c2) then
      n_bad_link = n_bad_link + 1
      cycle
    end if

    sx = Grid % sx(s);  sy = Grid % sy(s);  sz = Grid % sz(s)
    ss = sqrt(sx**2 + sy**2 + sz**2)
    if(ss .lt. TINY) then
      n_zero_area = n_zero_area + 1
      if(f_zero_area .eq. 0) f_zero_area = s
      cell_flag(c1) = 2.0
      if(c2 .gt. 0) cell_flag(c2) = 2.0
      cycle
    end if

    xc1 = Grid % xc(c1);  yc1 = Grid % yc(c1);  zc1 = Grid % zc(c1)

    !-----------------------!
    !   Boundary face       !
    !-----------------------!
    if(c2 .lt. 0) then
      n_bnd_faces = n_bnd_faces + 1

      b = Grid % region % at_cell(c2)
      if(b .lt. 1 .or. b .gt. Grid % n_bnd_regions) then
        n_bad_region = n_bad_region + 1
        cycle
      end if

      xd = Grid % xf(s) - xc1
      yd = Grid % yf(s) - yc1
      zd = Grid % zf(s) - zc1
      dd = sqrt(xd**2 + yd**2 + zd**2)

      ! Normal distance from cell center to the boundary face
      dn = abs(sx*xd + sy*yd + sz*zd) / ss
      b_count(b) = b_count(b) + 1
      b_area (b) = b_area (b) + ss
      b_dmin (b) = min(b_dmin(b), dn)
      b_dmax (b) = max(b_dmax(b), dn)
      b_dsum (b) = b_dsum (b) + dn

    !-----------------------!
    !   Inside face         !
    !-----------------------!
    else
      n_int_faces = n_int_faces + 1

      ! Volume ratio (valid for periodic faces too)
      vol1 = Grid % vol(c1)
      vol2 = Grid % vol(c2)
      if(vol1 .gt. 0.0 .and. vol2 .gt. 0.0) then
        ratio = max(vol1 / vol2, vol2 / vol1)
        n_vr   = n_vr + 1
        vr_sum = vr_sum + ratio
        if(ratio .gt. vr_max) then
          vr_max = ratio
          vr_max_at = (/Grid % xf(s), Grid % yf(s), Grid % zf(s)/)
        end if
        if(ratio .gt. VOL_RAT_WARN) then
          w_vr = w_vr + 1
          cell_flag(c1) = max(cell_flag(c1), 1.0)
          cell_flag(c2) = max(cell_flag(c2), 1.0)
        end if
      end if

      if(Grid % faces_s(s) .ne. 0) then
        n_per_faces = n_per_faces + 1
        cycle              ! periodic face, geometry not assessed
      end if

      xd = Grid % xc(c2) - xc1
      yd = Grid % yc(c2) - yc1
      zd = Grid % zc(c2) - zc1
      dd = sqrt(xd**2 + yd**2 + zd**2)
    end if

    !-------------------------------------------------!
    !   Non-orthogonality, for boundary and inside    !
    !-------------------------------------------------!
    if(dd .lt. TINY) then
      n_inverted = n_inverted + 1
      cycle
    end if
    ds = sx*xd + sy*yd + sz*zd
    if(ds .le. 0.0) then
      n_inverted = n_inverted + 1         ! cells' connection points backwards
      cell_flag(c1) = 2.0
      if(c2 .gt. 0) cell_flag(c2) = 2.0
      cycle
    end if
    cs    = min(1.0, ds / (ss * dd))
    angle = acos(cs) * 180.0 / PI
    n_no   = n_no + 1
    no_sum = no_sum + angle
    no_min = min(no_min, angle)
    if(angle .gt. no_max) then
      no_max = angle
      no_max_at = (/Grid % xf(s), Grid % yf(s), Grid % zf(s)/)
    end if
    i = 1 + count(angle .gt. bin_edge_no)
    hist_no(i) = hist_no(i) + 1
    if(angle .gt. NON_ORTH_WARN) then
      w_no = w_no + 1
      cell_flag(c1) = max(cell_flag(c1), 1.0)
      if(c2 .gt. 0) cell_flag(c2) = max(cell_flag(c2), 1.0)
    end if

    !----------------------------------!
    !   Skewness (inside faces only)   !
    !----------------------------------!
    if(c2 .gt. 0) then
      t   = (sx * (Grid % xf(s) - xc1)     &
           + sy * (Grid % yf(s) - yc1)     &
           + sz * (Grid % zf(s) - zc1)) / ds
      x_i = xc1 + t * xd
      y_i = yc1 + t * yd
      z_i = zc1 + t * zd
      skew = sqrt(  (Grid % xf(s) - x_i)**2    &
                  + (Grid % yf(s) - y_i)**2    &
                  + (Grid % zf(s) - z_i)**2) / dd
      n_sk   = n_sk + 1
      sk_sum = sk_sum + skew
      sk_min = min(sk_min, skew)
      if(skew .gt. sk_max) then
        sk_max = skew
        sk_max_at = (/Grid % xf(s), Grid % yf(s), Grid % zf(s)/)
      end if
      i = 1 + count(skew .gt. bin_edge_sk)
      hist_sk(i) = hist_sk(i) + 1
      if(skew .gt. SKEW_WARN) then
        w_sk = w_sk + 1
        cell_flag(c1) = max(cell_flag(c1), 1.0)
        cell_flag(c2) = max(cell_flag(c2), 1.0)
      end if
    end if
  end do

  !------------------------------------------------!
  !                                                !
  !   Pass through cells: aspect ratio, regularity !
  !                                                !
  !------------------------------------------------!
  n_vol_bad = 0;  n_open = 0;  n_inward = 0;  n_few_faces = 0;  n_shape = 0
  n_bad_node = 0
  f_vol_bad = 0;  f_open = 0;  f_inward = 0;  f_few = 0;  f_shape = 0
  n_ar = 0;  ar_max = 1.0;  ar_sum = 0.0;  w_ar = 0;  ar_max_at = 0.0
  min_vol = +HUGE;  max_vol = -HUGE;  tot_vol = 0.0;  min_vol_at = 0.0

  allocate(node_used(Grid % n_nodes));  node_used(:) = .false.

  do c = 1, nc

    ! Volumes
    tot_vol = tot_vol + Grid % vol(c)
    max_vol = max(max_vol, Grid % vol(c))
    if(Grid % vol(c) .lt. min_vol) then
      min_vol = Grid % vol(c)
      min_vol_at = (/Grid % xc(c), Grid % yc(c), Grid % zc(c)/)
    end if
    if(Grid % vol(c) .le. 0.0) then
      n_vol_bad = n_vol_bad + 1
      if(f_vol_bad .eq. 0) f_vol_bad = c
      cell_flag(c) = 2.0
    end if

    ! Nodes of the cell
    n = abs(Grid % cells_n_nodes(c))
    do i = 1, n
      if(Grid % cells_n(i,c) .lt. 1 .or.  &
         Grid % cells_n(i,c) .gt. Grid % n_nodes) then
        n_bad_node = n_bad_node + 1
      else
        node_used(Grid % cells_n(i,c)) = .true.
      end if
    end do

    ! Number of faces against the shape
    n = Grid % cells_n_faces(c)
    if(n .lt. 4) then
      n_few_faces = n_few_faces + 1
      if(f_few .eq. 0) f_few = c
      cell_flag(c) = 2.0
      cycle
    end if
    bad_cell = .false.
    select case(Grid % cells_n_nodes(c))
      case(4);  bad_cell = (n .ne. 4)   ! tetrahedron
      case(5);  bad_cell = (n .ne. 5)   ! pyramid
      case(6);  bad_cell = (n .ne. 5)   ! wedge
      case(8);  bad_cell = (n .ne. 6)   ! hexahedron
    end select
    if(bad_cell) then
      n_shape = n_shape + 1
      if(f_shape .eq. 0) f_shape = c
      cell_flag(c) = max(cell_flag(c), 1.0)
    end if

    ! Faces of the cell: closure, orientation and aspect ratio
    sum_x = 0.0;  sum_y = 0.0;  sum_z = 0.0;  sum_s = 0.0
    d_min = +HUGE;  d_max = -HUGE
    bad_cell = .false.
    do i_fac = 1, n
      s = Grid % cells_f(i_fac, c)

      if    (Grid % faces_c(1,s) .eq. c) then
        s_sign = +1.0
      else if(Grid % faces_c(2,s) .eq. c) then
        s_sign = -1.0
      else
        n_bad_link = n_bad_link + 1
        cycle
      end if

      sum_x = sum_x + s_sign * Grid % sx(s)
      sum_y = sum_y + s_sign * Grid % sy(s)
      sum_z = sum_z + s_sign * Grid % sz(s)
      sum_s = sum_s + sqrt(  Grid % sx(s)**2 + Grid % sy(s)**2    &
                           + Grid % sz(s)**2)

      ! Face which physically encloses the cell (shadow, if periodic)
      sf = s
      if(Grid % faces_s(s) .ne. 0) then
        if(.not. Grid % Is_Face_In_Cell(s, c)) sf = Grid % faces_s(s)
      end if

      xd = Grid % xf(sf) - Grid % xc(c)
      yd = Grid % yf(sf) - Grid % yc(c)
      zd = Grid % zf(sf) - Grid % zc(c)
      dist = sqrt(xd**2 + yd**2 + zd**2)
      d_min = min(d_min, dist)
      d_max = max(d_max, dist)

      ! Outward orientation (not meaningful for concave cells)
      if(.not. Grid % concave(c)) then
        if(s_sign * (  Grid % sx(s) * xd + Grid % sy(s) * yd  &
                     + Grid % sz(s) * zd) .le. 0.0) bad_cell = .true.
      end if
    end do

    if(bad_cell) then
      n_inward = n_inward + 1
      if(f_inward .eq. 0) f_inward = c
      cell_flag(c) = 2.0
    end if

    if(sum_s .gt. 0.0) then
      if(sqrt(sum_x**2 + sum_y**2 + sum_z**2) / sum_s .gt. CLOSURE_TOL) then
        n_open = n_open + 1
        if(f_open .eq. 0) f_open = c
        cell_flag(c) = 2.0
      end if
    end if

    if(d_min .gt. TINY) then
      aspect = d_max / d_min
      n_ar   = n_ar + 1
      ar_sum = ar_sum + aspect
      if(aspect .gt. ar_max) then
        ar_max = aspect
        ar_max_at = (/Grid % xc(c), Grid % yc(c), Grid % zc(c)/)
      end if
      if(aspect .gt. ASPECT_WARN) then
        w_ar = w_ar + 1
        cell_flag(c) = max(cell_flag(c), 1.0)
      end if
    end if
  end do

  n_unused = count(.not. node_used)

  !----------------------------------------------------------!
  !   Connectivity: count domains which are not connected    !
  !----------------------------------------------------------!
  allocate(visited(nc));  visited(:) = .false.
  allocate(queue  (nc))
  n_comp   = 0
  min_comp = nc
  n_seen   = 0
  do i = 1, nc
    if(visited(i)) cycle
    n_comp = n_comp + 1
    head = 1;  tail = 1
    queue(1)   = i
    visited(i) = .true.
    n_seen = 0
    do while(head .le. tail)
      c = queue(head);  head = head + 1
      n_seen = n_seen + 1
      do i_fac = 1, Grid % cells_n_faces(c)
        s  = Grid % cells_f(i_fac, c)
        c2 = Grid % faces_c(1,s) + Grid % faces_c(2,s) - c
        if(c2 .ge. 1 .and. c2 .le. nc) then
          if(.not. visited(c2)) then
            visited(c2) = .true.
            tail = tail + 1
            queue(tail) = c2
          end if
        end if
      end do
    end do
    min_comp = min(min_comp, n_seen)
  end do

  !-------------------!
  !   Wall distance   !
  !-------------------!
  n_wd = 0;  n_wd_unset = 0;  wd_min = +HUGE;  wd_min_at = 0.0
  do c = 1, nc
    if(Grid % wall_dist(c) .gt. 0.0 .and.  &
       Grid % wall_dist(c) .lt. 0.5 * HUGE) then
      n_wd = n_wd + 1
      if(Grid % wall_dist(c) .lt. wd_min) then
        wd_min = Grid % wall_dist(c)
        wd_min_at = (/Grid % xc(c), Grid % yc(c), Grid % zc(c)/)
      end if
    else
      n_wd_unset = n_wd_unset + 1
    end if
  end do

  !--------------------------------------------!
  !                                            !
  !   Collect warnings and errors (messages)   !
  !                                            !
  !--------------------------------------------!
  if(n_vol_bad .gt. 0) call Add(.true., n_vol_bad,                            &
    'cells with zero or negative volume; first at cell', f_vol_bad)
  if(n_open .gt. 0) call Add(.true., n_open,                                  &
    'cells whose faces do not close; first at cell', f_open)
  if(n_inward .gt. 0) call Add(.true., n_inward,                              &
    'cells with inward-pointing faces; first at cell', f_inward)
  if(n_few_faces .gt. 0) call Add(.true., n_few_faces,                        &
    'cells with less than 4 faces; first at cell', f_few)
  if(n_zero_area .gt. 0) call Add(.true., n_zero_area,                        &
    'faces with zero area; first at face', f_zero_area)
  if(n_bad_link .gt. 0) call Add(.true., n_bad_link,                          &
    'faces with wrong cell links', 0)
  if(n_bad_region .gt. 0) call Add(.true., n_bad_region,                      &
    'boundary faces without a valid region', 0)
  if(n_inverted .gt. 0) call Add(.true., n_inverted,                          &
    'faces whose cell connection points against the face normal', 0)
  if(n_bad_node .gt. 0) call Add(.true., n_bad_node,                          &
    'invalid node indices in cells', 0)
  if(no_max .gt. NON_ORTH_BAD) call Add(.true., w_no,                         &
    'faces with non-orthogonality above 85 degrees', 0)
  if(sk_max .gt. SKEW_BAD) call Add(.true., w_sk,                             &
    'faces with skewness above 0.95', 0)

  if(n_shape .gt. 0) call Add(.false., n_shape,                               &
    'cells whose number of faces does not match the shape '  //               &
    '(hanging nodes?); first at cell', f_shape)
  if(no_max .gt. NON_ORTH_WARN .and. no_max .le. NON_ORTH_BAD)                &
    call Add(.false., w_no, 'faces with non-orthogonality above 70 degrees',0)
  if(sk_max .gt. SKEW_WARN .and. sk_max .le. SKEW_BAD)                        &
    call Add(.false., w_sk, 'faces with skewness above 0.8', 0)
  if(w_ar .gt. 0) call Add(.false., w_ar,                                     &
    'cells with aspect ratio above 1000', 0)
  if(w_vr .gt. 0) call Add(.false., w_vr,                                     &
    'faces between cells with volume ratio above 20', 0)
  if(n_unused .gt. 0) call Add(.false., n_unused,                             &
    'nodes which do not belong to any cell', 0)
  if(n_comp .gt. 1) call Add(.false., n_comp,                                 &
    'disconnected parts of the mesh (domains); smallest has cells', min_comp)
  if(n_wd_unset .gt. 0 .and. n_wd .gt. 0) call Add(.false., n_wd_unset,       &
    'cells without a valid wall distance', 0)

  !--------------------!
  !                    !
  !   Print report     !
  !                    !
  !--------------------!
  call Say('#=========================================================')
  call Say('#')
  call Say('# Mesh report for: '//trim(Grid % name))
  call Say('#')
  call Say('#---------------------------------------------------------')

  !---------------!
  !   Mesh size   !
  !---------------!
  call Grid % Bounding_Box(xmin, ymin, zmin, xmax, ymax, zmax)
  call Say('# Mesh size:')
  write(line, '(a,i12)') '#   Number of cells:            ', nc
  call Say(line)
  write(line, '(a,i12)') '#   Number of faces:            ', Grid % n_faces
  call Say(line)
  write(line, '(a,i12)') '#   Number of nodes:            ', Grid % n_nodes
  call Say(line)
  write(line, '(a,i12)') '#   Number of boundary faces:   ', n_bnd_faces
  call Say(line)
  write(line, '(a,i12)') '#   Number of periodic faces:   ', n_per_faces
  call Say(line)
  write(line, '(a,es12.4)') '#   Total volume:               ', tot_vol
  call Say(line)
  write(line, '(a,es12.4,a,3es10.2,a)') '#   Smallest cell volume:       ',    &
        min_vol, '  at (', min_vol_at, ')'
  call Say(line)
  write(line, '(a,es12.4)') '#   Largest cell volume:        ', max_vol
  call Say(line)
  write(line, '(a,3es11.3)') '#   Bounding box, minimum:     ',              &
        xmin, ymin, zmin
  call Say(line)
  write(line, '(a,3es11.3)') '#   Bounding box, maximum:     ',              &
        xmax, ymax, zmax
  call Say(line)
  call Say('#')

  !----------------------!
  !   Boundary regions   !
  !----------------------!
  write(line, '(a,i0)') '# Boundary regions: ', Grid % n_bnd_regions
  call Say(line)
  call Say('#   (distance is from the cell center to boundary face, '  //     &
           'measured along the face normal)')
  write(line, '(a4,1x,a20,a10,a12,3a11)') '# no', 'name', 'faces', 'area',    &
        'dist.min', 'dist.avg', 'dist.max'
  call Say(line)
  do b = 1, Grid % n_bnd_regions
    if(b_count(b) .gt. 0) then
      write(line, '(a2,i2,1x,a20,i10,es12.3,3es11.3)')                        &
            '# ', b, Grid % region % name(b)(1:20), b_count(b), b_area(b),    &
            b_dmin(b), b_dsum(b) / b_count(b), b_dmax(b)
    else
      write(line, '(a2,i2,1x,a20,i10,a)') '# ', b,                           &
            Grid % region % name(b)(1:20), 0, '  (no faces)'
    end if
    call Say(line)
  end do
  call Say('#')

  !-------------------!
  !   Wall distance   !
  !-------------------!
  call Say('# Wall distance:')
  if(n_wd .gt. 0) then
    write(line, '(a,es11.3)') '#   Minimum wall distance:      ', wd_min
    call Say(line)
    write(line, '(a,3es11.3,a)') '#   ... at cell center:        ',         &
          wd_min_at(1), wd_min_at(2), wd_min_at(3)
    call Say(line)
  else
    call Say('#   Not calculated (no wall selected or no walls found)')
  end if
  call Say('#')

  !----------------!
  !   Mesh quality !
  !----------------!
  call Say('# Mesh quality:')
  if(n_no .gt. 0) then
    write(line, '(a,3f9.2,a)') '#   Non-orthogonality (deg.)  min/avg/max:',  &
          no_min, no_sum / n_no, no_max, '   (limit 70)'
    call Say(line)
    write(line, '(a,3es10.2,a)') '#     worst face at:      ', no_max_at, ''
    call Say(line)
    write(line, '(a,6i9)') '#     histogram [<20,<40,<60,<70,<80,>=80]:',    &
          hist_no
    call Say(line)
  end if
  if(n_sk .gt. 0) then
    write(line, '(a,3f9.3,a)') '#   Skewness                  min/avg/max:',  &
          sk_min, sk_sum / n_sk, sk_max, '   (limit 0.8)'
    call Say(line)
    write(line, '(a,3es10.2,a)') '#     worst face at:      ', sk_max_at, ''
    call Say(line)
    write(line, '(a,6i9)') '#     histogram [<.1,<.25,<.5,<.8,<.95,>=.95]:', &
          hist_sk
    call Say(line)
  end if
  if(n_ar .gt. 0) then
    write(line, '(a,f9.2,a,f9.2)') '#   Cell aspect ratio         avg/max:  ',&
          ar_sum / n_ar, ' ', ar_max
    call Say(line)
    write(line, '(a,3es10.2,a)') '#     worst cell at:      ', ar_max_at, ''
    call Say(line)
  end if
  if(n_vr .gt. 0) then
    write(line, '(a,f9.2,a,f9.2)') '#   Volume ratio of neighbors avg/max:  ',&
          vr_sum / n_vr, ' ', vr_max
    call Say(line)
    write(line, '(a,3es10.2,a)') '#     worst face at:      ', vr_max_at, ''
    call Say(line)
  end if
  call Say('#')

  !----------------------!
  !   Regularity check   !
  !----------------------!
  call Say('# Cell regularity check:')
  call Say_Check('all cells have positive volume',            n_vol_bad)
  call Say_Check('all cells are closed by their faces',       n_open)
  call Say_Check('all faces point outwards of their cells',   n_inward)
  call Say_Check('all cells have at least 4 faces',           n_few_faces)
  call Say_Check('no degenerate (zero area) faces',           n_zero_area)
  call Say_Check('all faces correctly linked to their cells', n_bad_link)
  call Say_Check('all boundary faces belong to a region',     n_bad_region)
  call Say_Check('cells connect neighbors in face direction', n_inverted)
  call Say_Check('cell shapes consistent with face counts',   n_shape)
  call Say_Check('all nodes belong to some cell',             n_unused)
  call Say_Check('mesh is a single connected domain',         n_comp - 1)
  call Say('#')

  !------------------------------!
  !   Warnings and the verdict   !
  !------------------------------!
  if(n_err + n_warn .gt. 0) then
    call Say('# Findings:')
    do i = 1, n_err + n_warn
      call Say('# '//msg(i))
    end do
    call Say('#')
  end if

  if(n_err .gt. 0) then
    write(line, '(a,i0,a,i0,a)') '# VERDICT: mesh has ', n_err,               &
          ' error(s) and ', n_warn, ' warning(s).  Do not use it as it is!'
    call Say(line, RED)
  else if(n_warn .gt. 0) then
    write(line, '(a,i0,a)') '# VERDICT: mesh is regular, with ', n_warn,      &
          ' warning(s).  Please review them.'
    call Say(line, YELLOW)
  else
    call Say('# VERDICT: mesh is regular and no warnings were raised.', GREEN)
  end if
  call Say('#')
  call Say('#---------------------------------------------------------')

  close(fu)

  ! Save flagged cells for inspection in Paraview
  if(n_err + n_warn .gt. 0 .and. any(cell_flag .gt. 0.0)) then
    call Grid % Save_Debug_Vtu('mesh_report_flags',            &
                               scalar_cell = cell_flag,        &
                               scalar_name = 'cell_flag')
    print '(a)', ' # Cells with problems are flagged in file *mesh_report_'  //&
                 'flags* (cell_flag: 1 = warning, 2 = error)'
  end if

  call Profiler % Stop('Mesh_Report')

  contains

  !============================================================================!
  subroutine Say(text, color)
  !----------------------------------------------------------------------------!
  !   Prints a line on the screen and stores it in the report file             !
  !----------------------------------------------------------------------------!
    character(*), intent(in)           :: text
    character(*), intent(in), optional :: color
  !----------------------------------------------------------------------------!

    if(present(color)) then
      print '(a)', ' '//color//trim(text)//RESET
    else
      print '(a)', ' '//trim(text)
    end if
    write(fu, '(a)') trim(text)

  end subroutine

  !============================================================================!
  subroutine Say_Check(text, n_fail)
  !----------------------------------------------------------------------------!
  !   Prints the result of one regularity check                                !
  !----------------------------------------------------------------------------!
    character(*), intent(in) :: text
    integer,      intent(in) :: n_fail
  !----------------------------------------------------------------------------!

    if(n_fail .eq. 0) then
      call Say('#   [ OK ]   '//text, GREEN)
    else
      write(tmp, '(a,i0,a)') '   (', n_fail, ' found)'
      call Say('#   [FAILED] '//text//trim(tmp), RED)
    end if

  end subroutine

  !============================================================================!
  subroutine Add(is_error, n_items, text, at)
  !----------------------------------------------------------------------------!
  !   Stores a warning or an error to be listed at the end of the report       !
  !----------------------------------------------------------------------------!
    logical,      intent(in) :: is_error
    integer,      intent(in) :: n_items
    character(*), intent(in) :: text
    integer,      intent(in) :: at
  !----------------------------------------------------------------------------!
    character(256) :: buf
  !----------------------------------------------------------------------------!

    if(n_err + n_warn .ge. MAX_MSG) return

    if(at .gt. 0) then
      write(buf, '(a,i0,1x,a,1x,i0)') '#   ', n_items, text, at
    else
      write(buf, '(a,i0,1x,a)') '#   ', n_items, text
    end if

    if(is_error) then
      n_err = n_err + 1
      msg(n_err + n_warn) = 'ERROR:   '//trim(buf(5:))
    else
      n_warn = n_warn + 1
      msg(n_err + n_warn) = 'WARNING: '//trim(buf(5:))
    end if

  end subroutine

  end subroutine
