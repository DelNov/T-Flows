!==============================================================================!
  subroutine Calculate_Geometry(Generate, Grid, real_run)
!------------------------------------------------------------------------------!
!>  Calculates geometrical quantities (such as cell centers, face surface
!>  areas, ...) and handles periodicity of the grid.
!------------------------------------------------------------------------------!
!   This subroutine has a sibling in Convert_Mod, with the same name.  They    !
!   can never be quite the same, unfortunatelly, because the data they start   !
!   with is different.  One of the most distinct differences is the treatment  !
!   of periodicity.  Here, periodic faces are added to the existing (internal) !
!   ones, while in Convert_Mod existing faces are turned into periodic ones.   !
!                                                                              !
!   Functionality overview:                                                    !
!   * Geometry calculations: It computes geometric quantities such as cell     !
!     centers, face surface areas, face centers, and volumes of cells in grid. !
!   * Periodicity handling: The subroutine uniquely addresses the aspect of    !
!     periodicity in the grid. It adds periodic faces to the existing internal !
!     ones or adjusts existing faces to account for periodic boundaries.       !
!   * Boundary cell centers: It calculates the centers of boundary cells based !
!     on the geometry of the internal cells and their faces.                   !
!   * Cell volume calculation: It computes the volume of each cell             !
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
  class(Generate_Type) :: Generate  !! parent class
  type(Grid_Type)      :: Grid      !! grid being generated
  logical, intent(in)  :: real_run  !! logical variable, false for trial run,
                                    !! true for real run
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, m, s, n_per, nn, nf
  real    :: xs2, ys2, zs2, t, tot_surf
  integer :: fn(6,4)
!------------------------------------------------------------------------------!
  include 'Block_Numbering.h90'
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Generate)
!==============================================================================!

  call Profiler % Start('Calculate_Geometry')

  ! An error trap for c1 and c2
  do s = 1, Grid % n_faces
    if(Grid % faces_c(2,s) > 0) then
      if(Grid % faces_c(1,s) > Grid % faces_c(2,s)) then
        call Message % Error(60,                                        &
                 'This shoulnd''t have happened at real face! \n '  //  &
                 'This error is critical.  Exiting now.!')
      end if
    end if
  end do

  ! Copy face-node numbering for faces
  ! (Note that it is the same as for blocks here in Generator)
  fn = hex_block

  !--------------------------------------------!
  !   Define topology of all cells and faces   !
  !--------------------------------------------!
  do c = 1, Grid % n_cells
    Grid % cells_n_nodes(c) = 8
  end do
  do s = 1, Grid % n_faces
    Grid % faces_n_nodes(s) = 4
  end do

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

  !--------------------------------------!
  !   Calculate boundary cell centers    !
  !--------------------------------------!
  !   => depends on: xc,yc,zc,sx,sy,sz   !
  !   <= gives:      xc,yc,zc for c<0    !
  !--------------------------------------!
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
    end if
  end do ! through faces

  if(real_run) then

    !---------------------------------------------!
    !   Find the faces on the periodic boundary   !
    !---------------------------------------------!
    !   => depends on: xc,yc,zc,sx,sy,sz          !
    !   <= gives:      dx,dy,dz                   !
    !---------------------------------------------!

    ! Initialize variables for Grid periodicity
    n_per = 0
    Grid % per_x = 0.0
    Grid % per_y = 0.0
    Grid % per_z = 0.0

    do s = 1, Grid % n_faces

      ! Initialize
      Grid % dx(s) = 0.0
      Grid % dy(s) = 0.0
      Grid % dz(s) = 0.0

      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)
      if(c2 > 0) then

        ! Scalar product of the side with line c1 - c2 is a good ...
        ! ... criterion to check if face is on periodic boundary
        if( (Grid % sx(s) * (Grid % xc(c2) - Grid % xc(c1) ) +  &
             Grid % sy(s) * (Grid % yc(c2) - Grid % yc(c1) ) +  &
             Grid % sz(s) * (Grid % zc(c2) - Grid % zc(c1) )) < 0.0 ) then

          n_per = n_per + 1

          ! Find the coordinates of ...
          m = Grid % face_c_to_c(s,2)

          if(Grid % faces_n_nodes(s) .eq. 4) then

            ! Coordinates of the shadow face
            xs2=.25*(  Grid % xn(Grid % cells_n(fn(m,1), c2))  &
                     + Grid % xn(Grid % cells_n(fn(m,2), c2))  &
                     + Grid % xn(Grid % cells_n(fn(m,3), c2))  &
                     + Grid % xn(Grid % cells_n(fn(m,4), c2)))

            ys2=.25*(  Grid % yn(Grid % cells_n(fn(m,1), c2))  &
                     + Grid % yn(Grid % cells_n(fn(m,2), c2))  &
                     + Grid % yn(Grid % cells_n(fn(m,3), c2))  &
                     + Grid % yn(Grid % cells_n(fn(m,4), c2)))

            zs2=.25*(  Grid % zn(Grid % cells_n(fn(m,1), c2))  &
                     + Grid % zn(Grid % cells_n(fn(m,2), c2))  &
                     + Grid % zn(Grid % cells_n(fn(m,3), c2))  &
                     + Grid % zn(Grid % cells_n(fn(m,4), c2)))

            ! Store the shadow face nodes
            nn = Grid % faces_n_nodes(s)
            nf = Grid % n_faces + n_per   ! new face number
            Grid % faces_n_nodes(nf) = nn
            Grid % faces_n(1:nn, nf) = Grid % cells_n(fn(m,1:nn), c2)

          end if

          ! Work out the cell distance for periodic faces
          Grid % dx(s) = Grid % xf(s) - xs2  !------------------------!
          Grid % dy(s) = Grid % yf(s) - ys2  ! later: xc2 = xc2 + dx  !
          Grid % dz(s) = Grid % zf(s) - zs2  !------------------------!

          ! Store a few geometrical quantities for the new shadow face
          ! (Change the sign for surface areas (as the shadow face
          !  is oriented differently) but keep face distance the same.)
          Grid % xf(nf) =  xs2
          Grid % yf(nf) =  ys2
          Grid % zf(nf) =  zs2
          Grid % sx(nf) = -Grid % sx(s)
          Grid % sy(nf) = -Grid % sy(s)
          Grid % sz(nf) = -Grid % sz(s)
          Grid % dx(nf) =  Grid % dx(s)
          Grid % dy(nf) =  Grid % dy(s)
          Grid % dz(nf) =  Grid % dz(s)

          ! Store cells surrounding the face and each other's shadows
          Grid % faces_c(1, nf) = c1
          Grid % faces_c(2, nf) = c2
          Grid % faces_s(s)  = nf
          Grid % faces_s(nf) = s

          Grid % per_x = max(Grid % per_x, abs(Grid % dx(s)))
          Grid % per_y = max(Grid % per_y, abs(Grid % dy(s)))
          Grid % per_z = max(Grid % per_z, abs(Grid % dz(s)))

        end if !  s*(c2-c1) < 0.0
      end if  !  c2 > 0
    end do    !  faces

    ! More likely, one should do: Grid % n_shadows = Grid % n_shadows + n_per
    Grid % n_shadows = n_per

    print '(a38,i7)',   '# Number of periodic faces:          ', n_per
    print '(a38,f8.3)', '# Periodicity in x direction         ', Grid % per_x
    print '(a38,f8.3)', '# Periodicity in y direction         ', Grid % per_y
    print '(a38,f8.3)', '# Periodicity in z direction         ', Grid % per_z

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

    !------------------------------------!
    !   Calculate cell inertia tensors   !
    !------------------------------------!
    call Grid % Calculate_Cell_Inertia()

    !------------------------------------------------------------!
    !   Calculate the interpolation factors for the cell faces   !
    !------------------------------------------------------------!
    call Grid % Calculate_Face_Interpolation()
  end if  ! real_run

  call Profiler % Stop('Calculate_Geometry')

  end subroutine
