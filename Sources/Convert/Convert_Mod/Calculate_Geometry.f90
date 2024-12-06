!==============================================================================!
  subroutine Calculate_Geometry(Convert, Grid, ask, g)
!------------------------------------------------------------------------------!
!>  Calculates geometrical quantities (such as cell centers, face surface
!>  areas, ...) and handles periodicity of the grid.
!------------------------------------------------------------------------------!
!   This subroutine has a sibling in Generate_Mod, with the same name.  They   !
!   can never be quite the same, unfortunatelly, because the data they start   !
!   with is different. One of the most distinct differences is the treatment   !
!   of periodicity. Here, periodic faces are created from existing (internal)  !
!   ones, whereas in Generate_Mod, they are added to existing cells.           !
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initial checks and scaling: The subroutine begins with error checking    !
!     for cell face orientation and offers the option to scale the grid        !
!     geometry based on user input.                                            !
!   * Grid homogeneity: It assesses the homogeneity of the grid, potentially   !
!     enforcing uniformity for certain cases (important for polyhedral grids). !
!   * Cell center calculation: Calculates the centers of cells in the grid,    !
!     essential for subsequent geometrical computations.                       !
!   * Face surface area calculation: Computes the surface area of each face in !
!     the grid, a key step in volume and flux calculations.                    !
!   * Face center calculation: Determines the geometric center of each face in !
!     the grid.                                                                !
!   * Boundary cell center calculation: The subroutine calculates the centers  !
!     of boundary cells.                                                       !
!   * Periodicity handling: This is a significant part of the subroutine,      !
!     where it identifies and processes faces on periodic boundaries.          !
!     It includes:                                                             !
!     - Finding periodic direction vectors and rotating vectors to align with  !
!       the periodic faces.                                                    !
!     - Matching periodic faces and creating 'shadow' faces to handle periodic !
!       conditions.                                                            !
!     - Adjusting boundary cell indices and coordinates to accommodate         !
!       periodicity.                                                           !
!   * Face orientation correction: The subroutine checks and corrects the      !
!     orientation of faces in the grid to ensure consistency in calculation.   !
!   * Volume calculations: Computes the volume of each cell.                   !
!   * Cell inertia tensor calculation: Calculates the inertia tensor for each  !
!     cell, providing information about the cell's resistance to rotation.     !
!   * Face interpolation factors: Determines interpolation factors for the     !
!     cell faces, facilitating the estimation of properties at face centers    !
!     from cell center values.                                                 !
!------------------------------------------------------------------------------!
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
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
  logical, intent(in) :: ask      !! controls interactive or optional parts
  integer, intent(in) :: g        !! grid rank
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, n, n1, n2, s, b, i, j
  integer              :: c11, c12, c21, c22, s1, s2, cnt_bnd, cnt_per
  integer              :: reg_per, n_per, number_faces
  real                 :: xs2, ys2, zs2
  real                 :: t, tot_surf, dis, min_dis, max_dis
  real                 :: v(3), k(3), v_o(3), v_r(3), theta  ! for rotation
  real,    allocatable :: b_coor_1(:), b_coor_2(:), b_coor_3(:)
  integer, allocatable :: b_face(:)
  character(SL)        :: answer, bou_cen
  real                 :: factor, prod
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  call Profiler % Start('Calculate_Geometry')

  ! An error trap for c1 and c2
  do s = 1, Grid % n_faces
    if(Grid % faces_c(2,s) > 0) then
      if(Grid % faces_c(1,s) > Grid % faces_c(2,s)) then
        call Message % Error(60,                                       &
                 "This shoulnd't have happened at real face! \n "  //  &
                 "This error is critical.  Exiting now.!",             &
                 file=__FILE__, line=__LINE__)
      end if
    end if
  end do

  !-------------------------------!
  !   Scale geometry              !
  !-------------------------------!
  !   => depends on: xn, yn, zn   !
  !   <= gives:      xn, yn, zn   !
  !-------------------------------!
  if(ask) then
    print *, '#========================================='
    print *, '# Geometric extents:                 '
    print *, '#-----------------------------------------'
    print '(2(a,es10.3))', ' # X from: ', minval(Grid % xn(:)),  &
                           '  to: ',      maxval(Grid % xn(:))
    print '(2(a,es10.3))', ' # Y from: ', minval(Grid % yn(:)),  &
                           '  to: ',      maxval(Grid % yn(:))
    print '(2(a,es10.3))', ' # Z from: ', minval(Grid % zn(:)),  &
                           '  to: ',      maxval(Grid % zn(:))
    print *, '# Enter scaling factor for geometry: '
    print *, '# (or skip to keep as is): '
    print *, '#-----------------------------------------'
    call File % Read_Line(5)
    answer = Line % tokens(1)
    call String % To_Upper_Case(answer)

    if( answer .ne. 'SKIP' ) then
      read(Line % tokens(1), *) factor
      print '(a,es10.3)', ' # Scaling geometry by factor: ', factor
      Grid % xn(:) = Grid % xn(:) * factor
      Grid % yn(:) = Grid % yn(:) * factor
      Grid % zn(:) = Grid % zn(:) * factor
    end if
  end if

  ! Check homogeneity of the grid and ask if you want it uniform
  if(g.eq.1)                                                          &
    call Grid % Search_Coordinate_Clusters(nodal           = .true.,  &
                                           enforce_uniform = .false.)
  if(g.eq.2)                                                          &
    call Grid % Search_Coordinate_Clusters(nodal           = .true.,  &
                                           enforce_uniform = .true.)

  !-----------------------------------------!
  !   Calculate the cell centers            !
  !-----------------------------------------!
  !   => depends on: xn, yn, zn             !
  !   <= gives:      xc, yc, zc @ c > 0     !
  !-----------------------------------------!
  call Grid % Calculate_Cell_Centers()

  !----------------------------------!
  !   Calculate face surface areas   !
  !----------------------------------!
  !   => depends on: xn, yn, zn      !
  !   <= gives:      sx, sy, sz      !
  !----------------------------------!
  call Grid % Calculate_Face_Surfaces()

  !--------------------------------!
  !   Calculate the face centers   !
  !--------------------------------!
  !   => depends on: xn, yn, zn    !
  !   <= gives:      xf, yf, zf    !
  !--------------------------------!
  call Grid % Calculate_Face_Centers()

  !-------------------------------------------!
  !   Calculate boundary cell centers         !
  !-------------------------------------------!
  !   => depends on: xc, yc, zc, sx, sy, sz   !
  !   <= gives:      xc, yc, zc  for c<0      !
  !-------------------------------------------!
  if(ask) then
    print *, '#===================================='
    print *, '# Position the boundary cell centres:'
    print *, '#------------------------------------'
    print *, '# Type 1 for barycentric placement'
    print *, '# Type 2 for orthogonal placement'
    print *, '#------------------------------------'
    bou_cen = File % Single_Word_From_Keyboard()
  else
    bou_cen = '1'
  end if

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    tot_surf = sqrt(Grid % sx(s)*Grid % sx(s)    &
                  + Grid % sy(s)*Grid % sy(s)    &
                  + Grid % sz(s)*Grid % sz(s) )

    if(c2 < 0) then
      t = (   Grid % sx(s) * (Grid % xf(s) - Grid % xc(c1))        &
            + Grid % sy(s) * (Grid % yf(s) - Grid % yc(c1))        &
            + Grid % sz(s) * (Grid % zf(s) - Grid % zc(c1)) ) / tot_surf
      Grid % xc(c2) = Grid % xc(c1) + Grid % sx(s)*t / tot_surf
      Grid % yc(c2) = Grid % yc(c1) + Grid % sy(s)*t / tot_surf
      Grid % zc(c2) = Grid % zc(c1) + Grid % sz(s)*t / tot_surf
      if(bou_cen .eq. '1') then
        Grid % xc(c2) = Grid % xf(s)
        Grid % yc(c2) = Grid % yf(s)
        Grid % zc(c2) = Grid % zf(s)
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
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 > 0) then
      Grid % dx(c1) = max( Grid % dx(c1), abs( Grid % xc(c2) - Grid % xc(c1) ) )
      Grid % dy(c1) = max( Grid % dy(c1), abs( Grid % yc(c2) - Grid % yc(c1) ) )
      Grid % dz(c1) = max( Grid % dz(c1), abs( Grid % zc(c2) - Grid % zc(c1) ) )
      Grid % dx(c2) = max( Grid % dx(c2), abs( Grid % xc(c2) - Grid % xc(c1) ) )
      Grid % dy(c2) = max( Grid % dy(c2), abs( Grid % yc(c2) - Grid % yc(c1) ) )
      Grid % dz(c2) = max( Grid % dz(c2), abs( Grid % zc(c2) - Grid % zc(c1) ) )
    end if
  end do ! through faces

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      if( Math % Approx_Real(Grid % dx(c1), 0.0) )  &
        Grid % xc(c1) = 0.75 * Grid % xc(c1) + 0.25 * Grid % xc(c2)
      if( Math % Approx_Real(Grid % dy(c1), 0.0) )  &
        Grid % yc(c1) = 0.75 * Grid % yc(c1) + 0.25 * Grid % yc(c2)
      if( Math % Approx_Real(Grid % dz(c1), 0.0) )  &
        Grid % zc(c1) = 0.75 * Grid % zc(c1) + 0.25 * Grid % zc(c2)
    end if
  end do ! through faces

  ! Why are the following three lines needed?
  ! Because memory for dx, dy and dz was used in the previous step
  Grid % dx(:) = 0.0
  Grid % dy(:) = 0.0
  Grid % dz(:) = 0.0

  !-----------------------------------------------------!
  !   For polyhedral grids, faces can be messy.         !
  !   This is the first place where they can be fixed   !
  !-----------------------------------------------------!
  !   => depends on:                                    !
  !   <= gives:                                         !
  !-----------------------------------------------------!
  call Grid % Correct_Face_Surfaces()

  !--------------------------------------------!
  !   Find the faces on the periodic boundary  !
  !--------------------------------------------!
  !   => depends on: xc, yc, zc, sx, sy, sz    !
  !   <= gives:      dx, dy, dz                !
  !--------------------------------------------!
  allocate(b_coor_1(Grid % n_faces)); b_coor_1 = 0.0
  allocate(b_coor_2(Grid % n_faces)); b_coor_2 = 0.0
  allocate(b_coor_3(Grid % n_faces)); b_coor_3 = 0.0
  allocate(b_face  (Grid % n_faces)); b_face   = 0

  !--------------------------------------------------------!
  !                                                        !
  !   Phase I  ->  find the faces on periodic boundaries   !
  !                                                        !
  !--------------------------------------------------------!
  if(ask) then
    answer = ''
    do while(answer .ne. 'SKIP')

      call Grid % Print_Regions_List()
      n_per = 0
      print *, '#=============================================================='
      print *, '# Enter the ordinal number(s) of periodic-boundary condition(s)'
      print *, '# from the boundary condition list (see above)                 '
      print *, '# Type skip if there is none !                                 '
      print *, '#--------------------------------------------------------------'
      call File % Read_Line(5)
      answer = Line % tokens(1)
      call String % To_Upper_Case(answer)

      if( answer .eq. 'SKIP' ) then
        reg_per = 0
        exit
      end if

      read(Line % tokens(1), *) reg_per
      if( reg_per > Grid % n_bnd_regions ) then
        print *, '# Critical error: boundary condition ', reg_per,  &
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
      do s = 1, Grid % n_faces
        c2 = Grid % faces_c(2,s)
        if(c2 < 0) then
          if(Grid % region % at_cell(c2) .eq. reg_per) then
            cnt_per = cnt_per + 1

            ! This is a dot product of surface vector and vector 1.0, 1.0, 1,0
            if( Grid % sx(s) + Grid % sy(s) + Grid % sz(s) > 0.0 ) then
              v(1) = v(1) + Grid % sx(s)
              v(2) = v(2) + Grid % sy(s)
              v(3) = v(3) + Grid % sz(s)
            else
              v(1) = v(1) - Grid % sx(s)
              v(2) = v(2) - Grid % sy(s)
              v(3) = v(3) - Grid % sz(s)
            end if
          end if
        end if
      end do
      v(1:3) = v(1:3) / cnt_per;  v(1:3) = v(1:3) / norm2(v(1:3))
      k(1:3) = Math % Cross_Product(v(1:3), (/1.,0.,0./))
      theta  = acos(dot_product    (v(1:3), (/1.,0.,0./)))
      print '(a)',       ' #==================================================='
      print '(a,3f7.3)', ' # Periodic direction vector: ', v(1:3)
      print '(a,3f7.3)', ' # Rotational vector:       : ', k(1:3)
      print '(a,3f7.3)', ' # Rotational angle:        : ', theta * 57.2957795131
      print '(a)',       ' #---------------------------------------------------'

      !---------------------------------------------------!
      !   Fill up helping vectors with sorting criteria   !
      !---------------------------------------------------!
      cnt_per = 0
      do s = 1, Grid % n_faces
        c2 = Grid % faces_c(2,s)
        if(c2 < 0) then
          if(Grid % region % at_cell(c2) .eq. reg_per) then
            v_o(1) = Grid % xf(s)
            v_o(2) = Grid % yf(s)
            v_o(3) = Grid % zf(s)
            v_r(1:3) = Math % Rotate_Vector(v_o(1:3), k(1:3), theta)
            cnt_per = cnt_per + 1
            b_coor_1(cnt_per) = v_r(1)
            b_coor_2(cnt_per) = v_r(2)
            b_coor_3(cnt_per) = v_r(3)
            b_face(cnt_per)   = s
          end if
        end if
      end do

      !-------------------------------------------!
      !   Sort the faces at periodic boundaries   !
      !-------------------------------------------!
      call Sort % Three_Real_Carry_Int(b_coor_1(1:cnt_per),  &
                                       b_coor_2(1:cnt_per),  &
                                       b_coor_3(1:cnt_per),  &
                                       b_face(1:cnt_per))

      !---------------------------------------------!
      !   Match the periodic faces with shadows &   !
      !    fill up the Grid % faces_s structure     !
      !---------------------------------------------!

      ! Check if this simple sort worked out
      min_dis =  HUGE
      max_dis = -HUGE
      do s = 1, cnt_per / 2
        s1 = b_face(s)
        s2 = b_face(s + cnt_per / 2)
        dis = Math % Distance(Grid % xf(s1), Grid % yf(s1), Grid % zf(s1),  &
                              Grid % xf(s2), Grid % yf(s2), Grid % zf(s2))
        min_dis = min(min_dis, dis)
        max_dis = max(max_dis, dis)
      end do

      ! Simple sort didn't work out
      if( .not. Math % Approx_Real(min_dis, max_dis, tol=MICRO) ) then
        print *, '# Simple sort of periodic faces failed,'  //  &
                 ' trying a slower algorithm'
        do i = 1, cnt_per / 2
          s1 = b_face(i)
          do j = cnt_per / 2 + 1, cnt_per
            s2 = b_face(j)
            dis = Math % Distance(Grid % xf(s1), Grid % yf(s1), Grid % zf(s1), &
                                  Grid % xf(s2), Grid % yf(s2), Grid % zf(s2))
            if(Math % Approx_Real(min_dis, dis, tol=MICRO) ) then
              c11 = Grid % faces_c(1,s1)  ! cell 1 for face 1
              c21 = Grid % faces_c(2,s1)  ! cell 2 for cell 1
              c12 = Grid % faces_c(1,s2)  ! cell 1 for face 2
              c22 = Grid % faces_c(2,s2)  ! cell 2 for face 2
              Grid % faces_s(s1) = s2     ! store where it was coppied from ...
              Grid % faces_s(s2) = s1     ! ... and for the mirror face too
              Grid % faces_c(2,s1) = c12  ! inside cell on the other side
              Grid % faces_c(1,s2) = 0    ! c21; zero marks a shadow face, dirty
              Grid % faces_c(2,s2) = 0    ! c21; zero marks a shadow face, dirty
              goto 1
            end if
          end do  ! s2
1         continue
        end do    ! s1

      ! Simple sort worked out
      else
        print *, '# Simple sort for periodic faces worked out - good!'
        do s = 1, cnt_per / 2
          s1 = b_face(s)
          s2 = b_face(s + cnt_per / 2)
          c11 = Grid % faces_c(1,s1)  ! cell 1 for face 1
          c21 = Grid % faces_c(2,s1)  ! cell 2 for cell 1
          c12 = Grid % faces_c(1,s2)  ! cell 1 for face 2
          c22 = Grid % faces_c(2,s2)  ! cell 2 for face 2
          Grid % faces_s(s1) = s2     ! store where it was coppied from ...
          Grid % faces_s(s2) = s1     ! ... and for the mirror face too
          Grid % faces_c(2,s1) = c12  ! inside cell on the other side
          Grid % faces_c(1,s2) = 0    ! c21; this zero marks a shadow face, dirty
          Grid % faces_c(2,s2) = 0    ! c21; this zero marks a shadow face, dirty
        end do
      end if

      n_per = cnt_per / 2
      print *, '# Phase I: periodic cells: ', n_per

      !---------------------------------!
      !      Find periodic extents      !
      !   (This is actually obsolete)   !
      !---------------------------------!
      Grid % per_x = 0.0
      Grid % per_y = 0.0
      Grid % per_z = 0.0
      do s = 1, n_per
        s1 = b_face(s)
        s2 = Grid % faces_s(s1)
        Grid % per_x = max(Grid % per_x, abs(Grid % xf(s1) - Grid % xf(s2)))
        Grid % per_y = max(Grid % per_y, abs(Grid % yf(s1) - Grid % yf(s2)))
        Grid % per_z = max(Grid % per_z, abs(Grid % zf(s1) - Grid % zf(s2)))
      end do
      print '(a38,f8.3)', ' # Periodicity in x direction         ', Grid % per_x
      print '(a38,f8.3)', ' # Periodicity in y direction         ', Grid % per_y
      print '(a38,f8.3)', ' # Periodicity in z direction         ', Grid % per_z

      !-------------------------------------------------!
      !   Compress all boundary cells by removing all   !
      !   cells which were holding periodic condition   !
      !-------------------------------------------------!
      cnt_bnd = 0
      Grid % new_c = 0
      do c = -1, -Grid % n_bnd_cells, -1
        if(Grid % region % at_cell(c) .ne. reg_per) then
          cnt_bnd = cnt_bnd + 1
          Grid % new_c(c) = -cnt_bnd
        end if
      end do

      ! Compress coordinates
      do c = -1, -Grid % n_bnd_cells, -1
        if(Grid % new_c(c) .ne. 0) then
          Grid % xc(Grid % new_c(c)) = Grid % xc(c)
          Grid % yc(Grid % new_c(c)) = Grid % yc(c)
          Grid % zc(Grid % new_c(c)) = Grid % zc(c)
         Grid % region % at_cell(Grid % new_c(c)) = Grid % region % at_cell(c)
        end if
      end do

      ! Compress indices
      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        if(Grid % new_c(c2) .ne. 0) then
          Grid % faces_c(2,s) = Grid % new_c(c2)
        end if
      end do

      Grid % n_bnd_cells = cnt_bnd
      print *, '# Kept boundary cells: ', Grid % n_bnd_cells

      !------------------------------------------------------------------!
      !   Remove boundary condition with reg_per and compress the rest   !
      !------------------------------------------------------------------!
      if(reg_per < Grid % n_bnd_regions) then

        ! Set the regions of boundary selected to be periodic to zero
        do c = -1, -Grid % n_bnd_cells, -1
          if(Grid % region % at_cell(c) .eq. reg_per) then
            Grid % region % at_cell(c) = 0
          end if
        end do

        ! Shift the rest of the boundary cells
        do b = 1, Grid % n_bnd_regions - 1
          if(b .ge. reg_per) then

            ! Correct the names
            Grid % region % name(b) = Grid % region % name (b+1)

            ! Correct all boundary regions too
            do c = -1, -Grid % n_bnd_cells, -1
              if(Grid % region % at_cell(c) .eq. (b+1)) then
                Grid % region % at_cell(c) = b
              end if
            end do

          end if
        end do
      else
        Grid % region % name(Grid % n_bnd_regions) = ''
      end if
      Grid % n_bnd_regions = Grid % n_bnd_regions - 1

    end do  ! while answer .ne. 'SKIP'
  end if    ! ask

  !----------------------------------------------------!
  !                                                    !
  !   Phase II  ->  work out dx, dy and dz for faces   !
  !                                                    !
  !----------------------------------------------------!

  !----------------!
  !   Initialize   !
  !----------------!
  n_per = 0
  Grid % dx(:) = 0.0
  Grid % dy(:) = 0.0
  Grid % dz(:) = 0.0

  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then

      !-------------------------!
      !   Find periodic faces   !
      !-------------------------!
      if(Grid % faces_s(s) .ne. 0) then

        n_per = n_per + 1

        ! Find the coordinates of the shadow face
        xs2 = Grid % xf(Grid % faces_s(s))
        ys2 = Grid % yf(Grid % faces_s(s))
        zs2 = Grid % zf(Grid % faces_s(s))

        Grid % dx(s) = Grid % xf(s) - xs2  !-----------------------!
        Grid % dy(s) = Grid % yf(s) - ys2  ! later: xc2 = xc2 + dx !
        Grid % dz(s) = Grid % zf(s) - zs2  !-----------------------!

        Grid % dx(Grid % faces_s(s)) = Grid % dx(s)
        Grid % dy(Grid % faces_s(s)) = Grid % dy(s)
        Grid % dz(Grid % faces_s(s)) = Grid % dz(s)

      end if !  s*(c2-c1) < 0.0
    end if   !  c2 > 0
  end do     !  faces

  ! Should this maybe be:
  ! Grid % n_shadows = Grid % n_shadows + n_per ?
  ! Actually no, because n_per is being re-counted in the above loop.
  Grid % n_shadows = n_per

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

  ! At this point, Grid % n_faces includes shadow faces too
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    !-------------------------------------------------------------------------!
    !   Product of centres connection and surface normal should be positive   !
    !-------------------------------------------------------------------------!
    prod = Grid % sx(s) * (Grid % xc(c2) + Grid % dx(s) - Grid % xc(c1) )  &
         + Grid % sy(s) * (Grid % yc(c2) + Grid % dy(s) - Grid % yc(c1) )  &
         + Grid % sz(s) * (Grid % zc(c2) + Grid % dz(s) - Grid % zc(c1) )

    !----------------------------------------------------------!
    !   If it is not, change the orientations of the surface   !
    !   (Shadow faces, which are at this point included in     !
    !    the Grid % n_faces count, have c1 and c2 set to 0)    !
    !----------------------------------------------------------!
    if(prod < 0 .and. c1 > 0 .and. c2 > 0) then

      ! Increase the counters
      if(c2 .gt. 0) n1 = n1 + 1
      if(c2 .lt. 0) n2 = n2 + 1

      ! Reverse the order of face's nodes
      n = Grid % faces_n_nodes(s)  ! number of nodes in this face
      call Sort % Reverse_Order_Int(Grid % faces_n(1:n, s))

      ! Keep the first node first (important if it is concave)
      Grid % faces_n(1:n,s) = cshift(Grid % faces_n(1:n,s), -1)

      ! Change the orientation of calculated surface vector
      Grid % sx(s) = -Grid % sx(s)
      Grid % sy(s) = -Grid % sy(s)
      Grid % sz(s) = -Grid % sz(s)

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
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c1 > 0) then
      number_faces = number_faces  + 1
      Grid % new_f(s) = number_faces
    end if
  end do

  ! Assign numbers to shadow faces, these can only be iside
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c1 .eq. 0 .and. c2 .eq. 0) then  ! marked like that above -> dirty
      number_faces = number_faces  + 1
      Grid % new_f(s) = number_faces

      ! Restore cells surrounding it
      Grid % faces_c(1, s) = Grid % faces_c(1, Grid % faces_s(s))
      Grid % faces_c(2, s) = Grid % faces_c(2, Grid % faces_s(s))
    end if
  end do
  print '(a38,i9)', ' # Old number of faces:               ',  Grid % n_faces
  print '(a38,i9)', ' # New number of faces:               ',  number_faces

  !----------------------------------!
  !                                  !
  !   Phase IV  ->  sort the faces   !
  !                                  !
  !----------------------------------!
  call Grid % Sort_Faces_By_Index(Grid % new_f, Grid % n_faces)
  call Sort % Real_By_Index(Grid % n_faces, Grid % xf, Grid % new_f)
  call Sort % Real_By_Index(Grid % n_faces, Grid % yf, Grid % new_f)
  call Sort % Real_By_Index(Grid % n_faces, Grid % zf, Grid % new_f)
  call Sort % Real_By_Index(Grid % n_faces, Grid % sx, Grid % new_f)
  call Sort % Real_By_Index(Grid % n_faces, Grid % sy, Grid % new_f)
  call Sort % Real_By_Index(Grid % n_faces, Grid % sz, Grid % new_f)
  call Sort % Real_By_Index(Grid % n_faces, Grid % dx, Grid % new_f)
  call Sort % Real_By_Index(Grid % n_faces, Grid % dy, Grid % new_f)
  call Sort % Real_By_Index(Grid % n_faces, Grid % dz, Grid % new_f)

  ! Why not: Grid % n_faces = Grid % n_faces - Grid % n_shadows?
  Grid % n_faces = Grid % n_faces - n_per

  ! Final correction to shadow faces for Grid % faces_s and Grid faces_c
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % faces_s(s) > 0) then
      Grid % faces_s(s) = Grid % new_f(Grid % faces_s(s))
      Grid % faces_c(1, Grid % faces_s(s)) = Grid % faces_c(1, s)
      Grid % faces_c(2, Grid % faces_s(s)) = Grid % faces_c(2, s)
    end if
  end do

  !-----------------------------------!
  !   Check the periodic boundaries   !
  !-----------------------------------!
  max_dis = 0.0
  do s = 1, Grid % n_faces
    max_dis = max(max_dis, (  Grid % dx(s)*Grid % dx(s)  &
                            + Grid % dy(s)*Grid % dy(s)  &
                            + Grid % dz(s)*Grid % dz(s) ) )
  end do
  print '(a45,es12.5)', ' # Maximal distance of periodic boundary is: ',  &
                        sqrt(max_dis)

  !----------------------------------!
  !   Calculate the cell volumes     !
  !----------------------------------!
  !   => depends on: xc,yc,zc,       !
  !                  dx,dy,dz,       !
  !                  xsp, ysp, zsp   !
  !   <= gives:      vol             !
  !----------------------------------!
  call Grid % Calculate_Cell_Volumes()

  Grid % min_vol =  HUGE
  Grid % max_vol = -HUGE
  Grid % tot_vol = 0.0
  do c = 1, Grid % n_cells
    Grid % tot_vol = Grid % tot_vol + Grid % vol(c)
    Grid % min_vol = min(Grid % min_vol, Grid % vol(c))
    Grid % max_vol = max(Grid % max_vol, Grid % vol(c))
  end do
  print '(a45,es12.5)', ' # Minimal cell volume is:                   ',  &
        Grid % min_vol
  print '(a45,es12.5)', ' # Maximal cell volume is:                   ',  &
        Grid % max_vol
  print '(a45,es12.5)', ' # Total domain volume is:                   ',  &
        Grid % tot_vol
  print *, '# Cell volumes calculated !'

  if(Grid % min_vol < 0.0) then
    print *, '# Negative volume occured!'
    print *, '# Execution will halt now!'
    stop
  end if

  !------------------------------------!
  !   Calculate cell inertia tensors   !
  !------------------------------------!
  call Grid % Calculate_Cell_Inertia()

  !------------------------------------------------------------!
  !   Calculate the interpolation factors for the cell faces   !
  !------------------------------------------------------------!
  call Grid % Calculate_Face_Interpolation()

  call Profiler % Stop('Calculate_Geometry')

  end subroutine
