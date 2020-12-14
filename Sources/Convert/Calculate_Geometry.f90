!==============================================================================!
  subroutine Calculate_Geometry(grid)
!------------------------------------------------------------------------------!
!   Calculates geometrical quantities of the grid.                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, n, s, b
  integer              :: c11, c12, c21, c22, s1, s2, bou_cen, cnt_bnd, cnt_per
  integer              :: color_per, n_per, number_faces
  integer              :: n1, n2, i_nod, j_nod
  real                 :: xt(MAX_FACES_N_NODES),  &
                          yt(MAX_FACES_N_NODES),  &
                          zt(MAX_FACES_N_NODES)
  real                 :: xs2, ys2, zs2
  real                 :: x_cell_1, y_cell_1, z_cell_1, dv_1
  real                 :: x_cell_2, y_cell_2, z_cell_2, dv_2
  real                 :: t, sur_tot, max_dis
  real                 :: v(3), k(3), v_o(3), v_r(3), theta  ! for rotation
  real,    allocatable :: b_coor(:)
  integer, allocatable :: b_face(:)
  character(SL)        :: answer
  real                 :: big, small, factor, prod
!==============================================================================!
!                                                                              !
!                                n3                                            !
!                 +---------------!---------------+                            !
!                /|              /|              /|                            !
!               / |             / |             / |                            !
!              /  |          n2/  |            /  |                            !
!             +---------------!---------------+   |                            !
!             |   |           |   |           |   |                            !
!             |   |     o---- | s-------o     |   |                            !
!             |   +---- c1 ---|   !---- c2 ---|   +                            !
!             |  /            |  /n4          |  /                             !
!             | /             | /             | /                              !
!             |/              |/              |/                               !
!             +---------------!---------------+                                !
!                            n1                                                !
!                                                                              !
!   Notes:                                                                     !
!                                                                              !
!     ! face s is oriented from cell center c1 to cell center c2               !
!     ! c2 is greater then c1 inside the domain or smaller then 0              !
!       on the boundary                                                        !
!     ! nodes are denoted with n1 - n4                                         !
!                                                                              !
!            c3                                                                !
!             \  4-----3                                                       !
!              \/ \  . |                                                       !
!              /   \  +---c2                                                   !
!             /  .  \  |                                                       !
!            / .     \ |                                                       !
!           /.        \|                                                       !
!          1-----------2                                                       !
!                   |                                                          !
!                   c1                                                         !
!                                                                              !
!                                n3                                            !
!                 +---------------!-------+                                    !
!                /|            n2/|      /|                                    !
!               / |             !-------+ |                                    !
!              /  |            /|s|  c2 | |                                    !
!             +---------------+ | !-----| +                                    !
!             |   |           | |/n4    |/                                     !
!             |   |     c1    | !-------+                                      !
!             |   +-----------|n1 +                                            !
!             |  /            |  /                                             !
!             | /             | /                                              !
!             |/              |/                                               !
!             +---------------+                                                !
!                            n1                                                !
!                                                                              !
!------------------------------------------------------------------------------!
!   Generaly:                                                                  !
!                                                                              !
!   the equation of plane reads: A * x + B * y + C * z + D = 0                 !
!                                                                              !
!   and the equation of line:  x = x0 + t*rx                                   !
!                              y = y0 + t*ry                                   !
!                              z = z0 + t*rz                                   !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       !
!   In our case:                                                               !
!                                                                              !
!     line is a connection between the two cell centers:                       !
!                                                                              !
!     x = xc(c1) + t*(xc(c2)-xc(c1)) = xc(c1) + t*rx                           !
!     y = yc(c1) + t*(yc(c2)-yc(c1)) = yc(c1) + t*ry                           !
!     z = zc(c1) + t*(zc(c2)-zc(c1)) = zc(c1) + t*rz                           !
!                                                                              !
!                                                                              !
!     plane is a cell face:                                                    !
!                                                                              !
!     sx * x + sy * y + Sz * z = sx * xsp(s) + sy * ysp(s) + sz * zsp(s)       !
!                                                                              !
!     and the intersection is at:                                              !
!                                                                              !
!         sx*(xsp(s)-xc(c1)) + sy*(ysp(s)-yc(c1) + sz*(zsp(s)-zc(c1))          !
!     t = -----------------------------------------------------------          !
!                           rx*sx + ry*sy + rz*sz                              !
!                                                                              !
!==============================================================================!

  ! An error trap for c1 and c2
  do s = 1, grid % n_faces
    if(grid % faces_c(2,s) > 0) then
      if(grid % faces_c(1,s) > grid % faces_c(2,s)) then
        print *, '# TROUBLE: this shoulnd''t have happened at real face!'
        print *, '# This error is critical.  Exiting now.!'
        stop
      end if
    end if
  end do

  !-------------------------------!
  !   Scale geometry              !
  !-------------------------------!
  !   => depends on: xn, yn, zn   !
  !   <= gives:      xn, yn, zn   !
  !-------------------------------!
  print *, '#========================================='
  print *, '# Geometric extents:                 '
  print *, '#-----------------------------------------'
  print '(2(a,es10.3))', ' # X from: ', minval(grid % xn(:)),  &
                         '  to: ',      maxval(grid % xn(:))
  print '(2(a,es10.3))', ' # Y from: ', minval(grid % yn(:)),  &
                         '  to: ',      maxval(grid % yn(:))
  print '(2(a,es10.3))', ' # Z from: ', minval(grid % zn(:)),  &
                         '  to: ',      maxval(grid % zn(:))
  print *, '# Enter scaling factor for geometry: '
  print *, '# (or skip to keep as is): '
  print *, '#-----------------------------------------'
  call File_Mod_Read_Line(5)
  answer = line % tokens(1)
  call To_Upper_Case(answer)

  if( answer .ne. 'SKIP' ) then
    read(line % tokens(1), *) factor
    print '(a,es10.3)', ' # Scaling geometry by factor: ', factor
    grid % xn(:) = grid % xn(:) * factor
    grid % yn(:) = grid % yn(:) * factor
    grid % zn(:) = grid % zn(:) * factor
  end if

  ! Estimate big and small
  call Grid_Mod_Estimate_Big_And_Small(grid, big, small)

  !-----------------------------------------!
  !   Calculate the cell centers            !
  !-----------------------------------------!
  !   => depends on: xn, yn, zn             !
  !   <= gives:      xc, yc, zc @ c > 0     !
  !-----------------------------------------!
  call Grid_Mod_Calculate_Cell_Centers(grid)

  !----------------------------------!
  !   Calculate face surface areas   !
  !----------------------------------!
  !   => depends on: xn, yn, zn      !
  !   <= gives:      sx, sy, sz      !
  !----------------------------------!
  call Grid_Mod_Calculate_Face_Surfaces(grid)

  !--------------------------------!
  !   Calculate the face centers   !
  !--------------------------------!
  !   => depends on: xn, yn, zn    !
  !   <= gives:      xf, yf, zf    !
  !--------------------------------!
  call Grid_Mod_Calculate_Face_Centers(grid)

  !-------------------------------------------!
  !   Calculate boundary cell centers         !
  !-------------------------------------------!
  !   => depends on: xc, yc, zc, sx, sy, sz   !
  !   <= gives:      xc, yc, zc  for c<0      !
  !-------------------------------------------!
  print *, '#===================================='
  print *, '# Position the boundary cell centres:'
  print *, '#------------------------------------'
  print *, '# Type 1 for barycentric placement'
  print *, '# Type 2 for orthogonal placement'
  print *, '#------------------------------------'
  read(*,*) bou_cen

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    sur_tot = sqrt(  grid % sx(s)*grid % sx(s)  &
                   + grid % sy(s)*grid % sy(s)  &
                   + grid % sz(s)*grid % sz(s) )

    if(c2 < 0) then
      t = (   grid % sx(s)*(grid % xf(s) - grid % xc(c1))        &
            + grid % sy(s)*(grid % yf(s) - grid % yc(c1))        &
            + grid % sz(s)*(grid % zf(s) - grid % zc(c1)) ) / sur_tot
      grid % xc(c2) = grid % xc(c1) + grid % sx(s)*t / sur_tot
      grid % yc(c2) = grid % yc(c1) + grid % sy(s)*t / sur_tot
      grid % zc(c2) = grid % zc(c1) + grid % sz(s)*t / sur_tot
      if(bou_cen .eq. 1) then
        grid % xc(c2) = grid % xf(s)
        grid % yc(c2) = grid % yf(s)
        grid % zc(c2) = grid % zf(s)
      end if
    end if
  end do ! through faces

  !---------------------------------------------------------------!
  !   Move the centers of co-planar molecules towards the walls   !
  !---------------------------------------------------------------!
  !   => depends on: xc, yc, zc                                   !
  !   <= gives:      xc, yc, zc                                   !
  !   +  uses:       dx, dy, dz                                   !
  !---------------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 > 0) then
      grid % dx(c1) = max( grid % dx(c1), abs( grid % xc(c2) - grid % xc(c1) ) )
      grid % dy(c1) = max( grid % dy(c1), abs( grid % yc(c2) - grid % yc(c1) ) )
      grid % dz(c1) = max( grid % dz(c1), abs( grid % zc(c2) - grid % zc(c1) ) )
      grid % dx(c2) = max( grid % dx(c2), abs( grid % xc(c2) - grid % xc(c1) ) )
      grid % dy(c2) = max( grid % dy(c2), abs( grid % yc(c2) - grid % yc(c1) ) )
      grid % dz(c2) = max( grid % dz(c2), abs( grid % zc(c2) - grid % zc(c1) ) )
    end if
  end do ! through faces

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if( Math_Mod_Approx_Real(grid % dx(c1), 0.0, small) )  &
        grid % xc(c1) = 0.75 * grid % xc(c1) + 0.25 * grid % xc(c2)
      if( Math_Mod_Approx_Real(grid % dy(c1), 0.0, small) )  &
        grid % yc(c1) = 0.75 * grid % yc(c1) + 0.25 * grid % yc(c2)
      if( Math_Mod_Approx_Real(grid % dz(c1), 0.0, small) )  &
        grid % zc(c1) = 0.75 * grid % zc(c1) + 0.25 * grid % zc(c2)
    end if
  end do ! through faces

  ! Why are the following three lines needed?
  ! Because memory for dx, dy and dz was used in the previous step
  grid % dx(:) = 0.0
  grid % dy(:) = 0.0
  grid % dz(:) = 0.0

  !--------------------------------------------!
  !   Find the faces on the periodic boundary  !
  !--------------------------------------------!
  !   => depends on: xc, yc, zc, sx, sy, sz    !
  !   <= gives:      dx, dy, dz                !
  !--------------------------------------------!
  allocate(b_coor(grid % n_faces)); b_coor = 0.0
  allocate(b_face(grid % n_faces)); b_face = 0

  !--------------------------------------------------------!
  !                                                        !
  !   Phase I  ->  find the faces on periodic boundaries   !
  !                                                        !
  !--------------------------------------------------------!
  answer = ''
  do while(answer .ne. 'SKIP')

    call Grid_Mod_Print_Bnd_Cond_List(grid)
    n_per = 0
    print *, '#=============================================================='
    print *, '# Enter the ordinal number(s) of periodic-boundary condition(s)'
    print *, '# from the boundary condition list (see above)                 '
    print *, '# Type skip if there is none !                                 '
    print *, '#--------------------------------------------------------------'
    call File_Mod_Read_Line(5)
    answer = line % tokens(1)
    call To_Upper_Case(answer)

    if( answer .eq. 'SKIP' ) then
      color_per = 0
      exit
    end if

    read(line % tokens(1), *) color_per
    if( color_per > grid % n_bnd_cond ) then
      print *, '# Critical error: boundary condition ', color_per,  &
                 ' doesn''t exist!'
      print *, '# Exiting! '
      stop
    end if

    !-----------------------------!
    !   Find periodic direction   !
    !-----------------------------!
    cnt_per = 0
    v(1:3)  = 0.0
    ! Browse through all the faces at periodic bc
    ! and accumulate periodic direction vector
    do s = 1, grid % n_faces
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(grid % bnd_cond % color(c2) .eq. color_per) then
          cnt_per = cnt_per + 1

          ! This is a dot product of surface vector and vector 1.0, 1.0, 1,0
          if( grid % sx(s) + grid % sy(s) + grid % sz(s) > 0.0 ) then
            v(1:3) = v(1:3) + (/grid % sx(s), grid % sy(s), grid % sz(s)/)
          else
            v(1:3) = v(1:3) - (/grid % sx(s), grid % sy(s), grid % sz(s)/)
          end if
        end if
      end if
    end do
    v(1:3) = v(1:3) / cnt_per;  v(1:3) = v(1:3) / norm2(v(1:3))
    k(1:3) = Math_Mod_Cross_Product(v(1:3), (/1.,0.,0./))
    theta  = acos(dot_product      (v(1:3), (/1.,0.,0./)))
    print '(A)',       ' #====================================================='
    print '(A,3F7.3)', ' # Periodic direction vector: ', v(1:3)
    print '(A,3F7.3)', ' # Rotational vector:       : ', k(1:3)
    print '(A,3F7.3)', ' # Rotational angle:        : ', theta * 57.2957795131
    print '(A)',       ' #-----------------------------------------------------'

    !---------------------------------------------------!
    !   Fill up helping vectors with sorting criteria   !
    !---------------------------------------------------!
    cnt_per = 0
    do s = 1, grid % n_faces
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(grid % bnd_cond % color(c2) .eq. color_per) then
          v_o(1) = grid % xf(s)
          v_o(2) = grid % yf(s)
          v_o(3) = grid % zf(s)
          v_r(1:3) = Math_Mod_Rotate_Vector(v_o(1:3), k(1:3), theta)
          cnt_per = cnt_per + 1
          b_coor(cnt_per) = v_r(1)*big**2 + v_r(2)*big + v_r(3)
          b_face(cnt_per) = s
        end if
      end if
    end do

    !-------------------------------------------!
    !   Sort the faces at periodic boundaries   !
    !-------------------------------------------!
    call Sort_Mod_Real_Carry_Int(b_coor(1:cnt_per), b_face(1:cnt_per))

    !---------------------------------------------!
    !   Match the periodic faces with shadows &   !
    !    fill up the grid % faces_s structure     !
    !---------------------------------------------!
    do s = 1, cnt_per / 2
      s1 = b_face(s)
      s2 = b_face(s + cnt_per / 2)
      c11 = grid % faces_c(1,s1)  ! cell 1 for face 1
      c21 = grid % faces_c(2,s1)  ! cell 2 for cell 1
      c12 = grid % faces_c(1,s2)  ! cell 1 for face 2
      c22 = grid % faces_c(2,s2)  ! cell 2 for face 2
      grid % faces_s(s1) = s2     ! store where it was coppied from ...
      grid % faces_s(s2) = s1     ! ... and for the mirror face too
      grid % faces_c(2,s1) = c12  ! inside cell on the other side of periodicity
      grid % faces_c(1,s2) = 0    ! c21; this zero marks a shadow face -> dirty
      grid % faces_c(2,s2) = 0    ! c21; this zero marks a shadow face -> dirty
    end do

    n_per = cnt_per / 2
    print *, '# Phase I: periodic cells: ', n_per

    !---------------------------------!
    !      Find periodic extents      !
    !   (This is actually obsolete)   !
    !---------------------------------!
    grid % per_x = 0.0
    grid % per_y = 0.0
    grid % per_z = 0.0
    do s = 1, n_per
      s1 = b_face(s)
      s2 = b_face(s + n_per)
      grid % per_x = max(grid % per_x, abs(grid % xf(s1) - grid % xf(s2)))
      grid % per_y = max(grid % per_y, abs(grid % yf(s1) - grid % yf(s2)))
      grid % per_z = max(grid % per_z, abs(grid % zf(s1) - grid % zf(s2)))
    end do
    print '(a38,f8.3)', ' # Periodicity in x direction         ', grid % per_x
    print '(a38,f8.3)', ' # Periodicity in y direction         ', grid % per_y
    print '(a38,f8.3)', ' # Periodicity in z direction         ', grid % per_z

    !-------------------------------------------------!
    !   Compress all boundary cells by removing all   !
    !   cells which were holding periodic condition   !
    !-------------------------------------------------!
    cnt_bnd = 0
    grid % new_c = 0
    do c = -1, -grid % n_bnd_cells, -1
      if(grid % bnd_cond % color(c) .ne. color_per) then
        cnt_bnd = cnt_bnd + 1
        grid % new_c(c) = -cnt_bnd
      end if
    end do

    ! Compress coordinates
    do c = -1, -grid % n_bnd_cells, -1
      if(grid % new_c(c) .ne. 0) then
        grid % xc(grid % new_c(c)) = grid % xc(c)
        grid % yc(grid % new_c(c)) = grid % yc(c)
        grid % zc(grid % new_c(c)) = grid % zc(c)
       grid % bnd_cond % color(grid % new_c(c)) = grid % bnd_cond % color(c)
      end if
    end do

    ! Compress indices
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(grid % new_c(c2) .ne. 0) then
        grid % faces_c(2,s) = grid % new_c(c2)
      end if
    end do

    grid % n_bnd_cells = cnt_bnd
    print *, '# Kept boundary cells: ', grid % n_bnd_cells

    !--------------------------------------------------------------------!
    !   Remove boundary condition with color_per and compress the rest   !
    !--------------------------------------------------------------------!
    if(color_per < grid % n_bnd_cond) then

      ! Set the color of boundary selected to be periodic to zero
      do c = -1, -grid % n_bnd_cells, -1
        if(grid % bnd_cond % color(c) .eq. color_per) then
          grid % bnd_cond % color(c) = 0
        end if
      end do

      ! Shift the rest of the boundary cells
      do b = 1, grid % n_bnd_cond - 1
        if(b .ge. color_per) then

          ! Correct the names
          grid % bnd_cond % name(b) = grid % bnd_cond % name (b+1)

          ! Correct all boundary colors too
          do c = -1, -grid % n_bnd_cells, -1
            if(grid % bnd_cond % color(c) .eq. (b+1)) then
              grid % bnd_cond % color(c) = b
            end if
          end do

        end if
      end do
    else
      grid % bnd_cond % name(grid % n_bnd_cond) = ''
    end if
    grid % n_bnd_cond = grid % n_bnd_cond - 1

  end do  ! while answer .ne. 'SKIP'

  !----------------------------------------------------!
  !                                                    !
  !   Phase II  ->  work out dx, dy and dz for faces   !
  !                                                    !
  !----------------------------------------------------!

  !----------------!
  !   Initialize   !
  !----------------!
  n_per = 0
  grid % dx(:) = 0.0
  grid % dy(:) = 0.0
  grid % dz(:) = 0.0

  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 > 0) then

      !-------------------------!
      !   Find periodic faces   !
      !-------------------------!
      if(grid % faces_s(s) .ne. 0) then

        n_per = n_per + 1

        ! Find the coordinates of the shadow face
        xs2 = grid % xf(grid % faces_s(s))
        ys2 = grid % yf(grid % faces_s(s))
        zs2 = grid % zf(grid % faces_s(s))

        grid % dx(s) = grid % xf(s) - xs2  !-----------------------!
        grid % dy(s) = grid % yf(s) - ys2  ! later: xc2 = xc2 + dx !
        grid % dz(s) = grid % zf(s) - zs2  !-----------------------!

        grid % dx(grid % faces_s(s)) = grid % dx(s)
        grid % dy(grid % faces_s(s)) = grid % dy(s)
        grid % dz(grid % faces_s(s)) = grid % dz(s)

      end if !  s*(c2-c1) < 0.0
    end if   !  c2 > 0
  end do     !  faces

  ! Should this maybe be:
  ! grid % n_shadows = grid % n_shadows + n_per ?
  ! Actually no, because n_per is being re-counted in the above loop.
  grid % n_shadows = n_per

  print '(a38,i9)',   ' # Phase II: number of shadow faces:  ', n_per

  !----------------------------------------------------!
  !                                                    !
  !   For some files, seems to be those generated in   !
  !   SnappyMesh, face orientation seems to be quite   !
  !   the opposite from what I expect in inner cells   !
  !                                                    !
  !----------------------------------------------------!
  n1 = 0
  n2 = 0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    !-------------------------------------------------------------------------!
    !   Product of centres connection and surface normal should be positive   !
    !-------------------------------------------------------------------------!
    prod = grid % sx(s) * (grid % xc(c2) + grid % dx(s) - grid % xc(c1) )  &
         + grid % sy(s) * (grid % yc(c2) + grid % dy(s) - grid % yc(c1) )  &
         + grid % sz(s) * (grid % zc(c2) + grid % dz(s) - grid % zc(c1) )

    !----------------------------------------------------------!
    !   If it is not, change the orientations of the surface   !
    !----------------------------------------------------------!
    if(prod < 0) then

      ! Increase the counters
      if(c2 > 0) n1 = n1 + 1
      if(c2 < 0) n2 = n2 + 1

      ! Reverse the order of face's nodes
      n = grid % faces_n_nodes(s)  ! number of nodes in this face
      call Sort_Mod_Reverse_Order_Int(grid % faces_n(1:n, s))

      ! Keep the first node first (important if it is concave)
      grid % faces_n(1:n,s) = cshift(grid % faces_n(1:n,s), -1)

      ! Change the orientation of calculated surface vector
      grid % sx(s) = -grid % sx(s)
      grid % sy(s) = -grid % sy(s)
      grid % sz(s) = -grid % sz(s)

    end if
  end do

  !-------------------------------------------------!
  !   Print some info on changed face orientation   !
  !-------------------------------------------------!
  if(n1 .gt. 0 .or. n2 .gt. 0) then
    print '(a)',      ' #======================================================'
    print '(a,i9,a)', ' # Changed orientation of', n1, ' inner faces,'
    print '(a,i9,a)', ' #                    and', n2, ' boundary faces.'
    print '(a)',      ' #------------------------------------------------------'
  else
    print '(a)', ' #========================================================='
    print '(a)', ' # Checked that all faces had good orientation (c1 -> c2). '
    print '(a)', ' #---------------------------------------------------------'
  end if

  !-------------------------------------------------------!
  !                                                       !
  !   Phase III  ->  find the new numbers of cell faces   !
  !                                                       !
  !-------------------------------------------------------!
  number_faces = 0

  ! Assign numbers to real cells, inside and boundary
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c1 > 0) then
      number_faces = number_faces  + 1
      grid % new_f(s) = number_faces
    end if
  end do

  ! Assign numbers to shadow faces, these can only be iside
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)
    if(c1 .eq. 0 .and. c2 .eq. 0) then  ! marked like that above -> dirty
      number_faces = number_faces  + 1
      grid % new_f(s) = number_faces

      ! Restore cells surrounding it
      grid % faces_c(1, s) = grid % faces_c(1, grid % faces_s(s))
      grid % faces_c(2, s) = grid % faces_c(2, grid % faces_s(s))
    end if
  end do
  print '(a38,i9)', ' # Old number of faces:               ',  grid % n_faces
  print '(a38,i9)', ' # New number of faces:               ',  number_faces

  !----------------------------------!
  !                                  !
  !   Phase IV  ->  sort the faces   !
  !                                  !
  !----------------------------------!
  call Grid_Mod_Sort_Faces_By_Index(grid, grid % new_f, grid % n_faces)
  call Sort_Mod_Real_By_Index(grid % n_faces, grid % xf, grid % new_f)
  call Sort_Mod_Real_By_Index(grid % n_faces, grid % yf, grid % new_f)
  call Sort_Mod_Real_By_Index(grid % n_faces, grid % zf, grid % new_f)
  call Sort_Mod_Real_By_Index(grid % n_faces, grid % sx, grid % new_f)
  call Sort_Mod_Real_By_Index(grid % n_faces, grid % sy, grid % new_f)
  call Sort_Mod_Real_By_Index(grid % n_faces, grid % sz, grid % new_f)
  call Sort_Mod_Real_By_Index(grid % n_faces, grid % dx, grid % new_f)
  call Sort_Mod_Real_By_Index(grid % n_faces, grid % dy, grid % new_f)
  call Sort_Mod_Real_By_Index(grid % n_faces, grid % dz, grid % new_f)

  ! Why not: grid % n_faces = grid % n_faces - grid % n_shadows?
  grid % n_faces = grid % n_faces - n_per

  ! Final correction to shadow faces for grid % faces_s and grid faces_c
  do s = 1, grid % n_faces + grid % n_shadows
    if(grid % faces_s(s) > 0) then
      grid % faces_s(s) = grid % new_f(grid % faces_s(s))
      grid % faces_c(1, grid % faces_s(s)) = grid % faces_c(1, s)
      grid % faces_c(2, grid % faces_s(s)) = grid % faces_c(2, s)
    end if
  end do

  !-----------------------------------!
  !   Check the periodic boundaries   !
  !-----------------------------------!
  max_dis = 0.0
  do s = 1, grid % n_faces
    max_dis = max(max_dis, (  grid % dx(s)*grid % dx(s)  &
                            + grid % dy(s)*grid % dy(s)  &
                            + grid % dz(s)*grid % dz(s) ) )
  end do
  print '(a45,e12.5)', ' # Maximal distance of periodic boundary is: ',  &
                       sqrt(max_dis)

  !----------------------------------!
  !   Calculate the cell volumes     !
  !----------------------------------!
  !   => depends on: xc,yc,zc,       !
  !                  dx,dy,dz,       !
  !                  xsp, ysp, zsp   !
  !   <= gives:      vol             !
  !----------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    do i_nod = 1, grid % faces_n_nodes(s)  ! for all face types
      xt(i_nod) = grid % xn(grid % faces_n(i_nod,s))
      yt(i_nod) = grid % yn(grid % faces_n(i_nod,s))
      zt(i_nod) = grid % zn(grid % faces_n(i_nod,s))
    end do

    ! First cell
    x_cell_1 = grid % xc(c1)
    y_cell_1 = grid % yc(c1)
    z_cell_1 = grid % zc(c1)
    do i_nod = 1, grid % faces_n_nodes(s)  ! for all face types
      j_nod = i_nod + 1;  if(j_nod > grid % faces_n_nodes(s)) j_nod = 1

      dv_1 = Math_Mod_Tet_Volume(grid % xf(s), grid % yf(s), grid % zf(s),  &
                                 xt(i_nod),    yt(i_nod),    zt(i_nod),     &
                                 xt(j_nod),    yt(j_nod),    zt(j_nod),     &
                                 x_cell_1,     y_cell_1,     z_cell_1)
      grid % vol(c1) = grid % vol(c1) + abs(dv_1)
    end do  ! i_nod

    ! Second cell
    if(c2 > 0) then
      x_cell_2 = grid % xc(c2) + grid % dx(s)
      y_cell_2 = grid % yc(c2) + grid % dy(s)
      z_cell_2 = grid % zc(c2) + grid % dz(s)

      do i_nod = 1, grid % faces_n_nodes(s)  ! for all face types
        j_nod = i_nod + 1;  if(j_nod > grid % faces_n_nodes(s)) j_nod = 1

        dv_2 = Math_Mod_Tet_Volume(grid % xf(s), grid % yf(s), grid % zf(s),  &
                                   xt(i_nod),    yt(i_nod),    zt(i_nod),     &
                                   xt(j_nod),    yt(j_nod),    zt(j_nod),     &
                                   x_cell_2,     y_cell_2,     z_cell_2)
        grid % vol(c2) = grid % vol(c2) + abs(dv_2)
      end do  ! i_nod
    end if

  end do

  grid % min_vol =  HUGE
  grid % max_vol = -HUGE
  grid % tot_vol = 0.0
  do c = 1, grid % n_cells
    grid % tot_vol = grid % tot_vol + grid % vol(c)
    grid % min_vol = min(grid % min_vol, grid % vol(c))
    grid % max_vol = max(grid % max_vol, grid % vol(c))
  end do
  print '(a45,es12.5)', ' # Minimal cell volume is:                   ',  &
        grid % min_vol
  print '(a45,es12.5)', ' # Maximal cell volume is:                   ',  &
        grid % max_vol
  print '(a45,es12.5)', ' # Total domain volume is:                   ',  &
        grid % tot_vol
  print *, '# Cell volumes calculated !'

  if(grid % min_vol < 0.0) then
    print *, '# Negative volume occured!'
    print *, '# Execution will halt now!'
    stop
  end if

  deallocate(b_coor)
  deallocate(b_face)

  !------------------------------------------------------------!
  !   Calculate the interpolation factors for the cell faces   !
  !------------------------------------------------------------!
  call Grid_Mod_Calculate_Face_Interpolation(grid)

  print *, '# Interpolation factors calculated !'

  end subroutine
