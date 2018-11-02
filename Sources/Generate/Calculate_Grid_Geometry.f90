!==============================================================================!
  subroutine Calculate_Grid_Geometry(grid, rrun)
!------------------------------------------------------------------------------!
!   Calculates geometrical quantities of the grid.                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: HUGE, ONE_THIRD
  use Gen_Mod,   only: face_c_to_c
  use Grid_Mod,  only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)     :: grid
  logical, intent(in) :: rrun
!----------------------------------[Calling]-----------------------------------!
  real :: Tet_Volume
  real :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, m, n, s, n_per
  real    :: loc_x_node(4), loc_y_node(4), loc_z_node(4)
  real    :: x_cell_tmp, y_cell_tmp, z_cell_tmp    
  real    :: xs2, ys2, zs2
  real    :: dsc1, dsc2          !  for the interpolation factors
  real    :: t, tot_surf 
  real    :: xc1, yc1, zc1, xc2, yc2, zc2 
  real    :: x_min, x_max, y_min, y_max, z_min, z_max
  integer :: f4n(6,4)
  integer :: f3n(4,3)
!==============================================================================!
!
!                                n3 
!                 +---------------!---------------+
!                /|              /|              /|
!               / |             / |             / |
!              /  |          n2/  |            /  |
!             +---------------!---------------+   |
!             |   |           |   |           |   |
!             |   |     o---- | s-------o     |   |      
!             |   +---- c1 ---|   !---- c2 ---|   +       
!             |  /            |  /n4          |  /
!             | /             | /             | /
!             |/              |/              |/
!             +---------------!---------------+
!                            n1
!
!   Notes:
! 
!     ! side s is oriented from cell center c1 to cell center c2     
!     ! c2 is greater then c1 inside the domain or smaller then 0
!       on the boundary
!     ! nodes are denoted with n1 - n4
!
!            c3           
!             \  4-----3
!              \/ \  . |
!              /   \  +---c2
!             /  .  \  |
!            / .     \ |
!           /.        \|
!          1-----------2
!                   |
!                   c1
!
!                                n3 
!                 +---------------!-------+         
!                /|            n2/|      /|
!               / |             !-------+ |
!              /  |            /|s|  c2 | |
!             +---------------+ | !-----| +
!             |   |           | |/n4    |/
!             |   |     c1    | !-------+
!             |   +-----------|n1 +
!             |  /            |  /
!             | /             | /
!             |/              |/ 
!             +---------------+
!                            n1
! 
!------------------------------------------------------------------------------!
!   Generaly:
!
!   the equation of plane reads: A * x + B * y + C * z + D = 0
!
!   and the equation of line:  x = x0 + t*rx
!                              y = y0 + t*ry
!                              z = z0 + t*rz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   In our case:
!
!     line is a connection between the two cell centers:
!
!     x = xc(c1) + t*(xc(c2)-xc(c1)) = xc(c1) + t*rx
!     y = yc(c1) + t*(yc(c2)-yc(c1)) = yc(c1) + t*ry
!     z = zc(c1) + t*(zc(c2)-zc(c1)) = zc(c1) + t*rz
!    
!
!     plane is a cell face: 
!
!     Sx * x + Sy * y + Sz * z = Sx * xsp(s) + Sy * ysp(s) + Sz * zsp(s)
!  
!     and the intersection is at:
!  
!         Sx*(xsp(s)-xc(c1)) + Sy*(ysp(s)-yc(c1) + Sz*(zsp(s)-zc(c1)) 
!     t = -----------------------------------------------------------
!                           rx*Sx + ry*Sy + rz*Sz
!  
!------------------------------------------------------------------------------!
  data    f4n / 1, 1, 2, 4, 3, 5,     &
                2, 5, 6, 8, 7, 7,     &
                4, 6, 8, 7, 5, 8,     &
                3, 2, 4, 3, 1, 6  /

  data    f3n / 1,  1,  2,  3,        &
                2,  4,  4,  4,        &
                3,  2,  3,  1 /

    ! Without the following six lines, this procedure works for any grid
    do c = 1, grid % n_cells
      grid % cells_n_nodes(c)=8
    end do
    do s = 1, grid % n_faces
      grid % faces_n_nodes(s)=4 
    end do

    !-----------------------------------------!
    !   Calculate the cell centers            !
    !-----------------------------------------!
    !   => depends on: x_node,y_node,z_node   ! 
    !   <= gives:      xc,yc,zc c>0           !
    !-----------------------------------------!
    do c = 1, grid % n_cells
      grid % xc(c)=0.0
      grid % yc(c)=0.0
      grid % zc(c)=0.0
      do n = 1, grid % cells_n_nodes(c)
        grid % xc(c) = grid % xc(c) + grid % xn(grid % cells_n(n,c))  &
                     / (1.0*grid % cells_n_nodes(c))
        grid % yc(c) = grid % yc(c) + grid % yn(grid % cells_n(n,c))  &
                     / (1.0*grid % cells_n_nodes(c))
        grid % zc(c) = grid % zc(c) + grid % zn(grid % cells_n(n,c))  &
                     / (1.0*grid % cells_n_nodes(c))
      end do
    end do

    !-----------------------------------------!
    !   Calculate delta                       !
    !-----------------------------------------!
    !   => depends on: x_node,y_node,z_node   ! 
    !   <= gives:      delta                  !
    !-----------------------------------------!
    do c = 1, grid % n_cells
      grid % delta(c)=0.0
      x_min = +HUGE   
      y_min = +HUGE  
      z_min = +HUGE  
      x_max = -HUGE  
      y_max = -HUGE  
      z_max = -HUGE  
      do n = 1, grid % cells_n_nodes(c)
        x_min = min(x_min, grid % xn(grid % cells_n(n,c)))
        y_min = min(y_min, grid % yn(grid % cells_n(n,c)))
        z_min = min(z_min, grid % zn(grid % cells_n(n,c)))
        x_max = max(x_max, grid % xn(grid % cells_n(n,c)))
        y_max = max(y_max, grid % yn(grid % cells_n(n,c)))
        z_max = max(z_max, grid % zn(grid % cells_n(n,c)))
      end do
      grid % delta(c) = x_max-x_min
      grid % delta(c) = max(grid % delta(c), (y_max-y_min))
      grid % delta(c) = max(grid % delta(c), (z_max-z_min))
    end do

    !-----------------------------------------------------!
    !   Calculate:                                        ! 
    !      components of cell sides, cell side centers.   !
    !-----------------------------------------------------!
    !   => depends on: x_node,y_node,z_node               ! 
    !   <= gives:      Sx,Sy,Sz,xsp,yzp,zsp               !
    !-----------------------------------------------------!
    do s = 1, grid % n_faces
      do n = 1, grid % faces_n_nodes(s)  ! for quadrilateral an triangular faces
        loc_x_node(n) = grid % xn(grid % faces_n(n,s))
        loc_y_node(n) = grid % yn(grid % faces_n(n,s))
        loc_z_node(n) = grid % zn(grid % faces_n(n,s))
      end do                       

      ! Cell side components
      if( grid % faces_n_nodes(s) .eq. 4 ) then
        grid % sx(s)= 0.5 * ((loc_y_node(2)-loc_y_node(1))  &
                           * (loc_z_node(2)+loc_z_node(1))  &
                           + (loc_y_node(3)-loc_y_node(2))  &
                           * (loc_z_node(2)+loc_z_node(3))  &
                           + (loc_y_node(4)-loc_y_node(3))  &
                           * (loc_z_node(3)+loc_z_node(4))  &
                           + (loc_y_node(1)-loc_y_node(4))  &
                           * (loc_z_node(4)+loc_z_node(1)) )
        grid % sy(s)= 0.5 * ((loc_z_node(2)-loc_z_node(1))  &
                           * (loc_x_node(2)+loc_x_node(1))  &
                           + (loc_z_node(3)-loc_z_node(2))  &
                           * (loc_x_node(2)+loc_x_node(3))  &
                           + (loc_z_node(4)-loc_z_node(3))  &
                           * (loc_x_node(3)+loc_x_node(4))  &
                           + (loc_z_node(1)-loc_z_node(4))  &
                           * (loc_x_node(4)+loc_x_node(1)) )
        grid % sz(s)= 0.5 * ((loc_x_node(2)-loc_x_node(1))  &
                           * (loc_y_node(2)+loc_y_node(1))  &
                           + (loc_x_node(3)-loc_x_node(2))  &
                           * (loc_y_node(2)+loc_y_node(3))  &
                           + (loc_x_node(4)-loc_x_node(3))  &
                           * (loc_y_node(3)+loc_y_node(4))  &
                           + (loc_x_node(1)-loc_x_node(4))  &
                           * (loc_y_node(4)+loc_y_node(1)) )
      else if( grid % faces_n_nodes(s) .eq. 3 ) then 
        grid % sx(s)= 0.5 * ((loc_y_node(2)-loc_y_node(1))  &
                           * (loc_z_node(2)+loc_z_node(1))  &
                           + (loc_y_node(3)-loc_y_node(2))  &
                           * (loc_z_node(2)+loc_z_node(3))  &
                           + (loc_y_node(1)-loc_y_node(3))  &
                           * (loc_z_node(3)+loc_z_node(1)) )
        grid % sy(s)= 0.5 * ((loc_z_node(2)-loc_z_node(1))  &
                           * (loc_x_node(2)+loc_x_node(1))  &
                           + (loc_z_node(3)-loc_z_node(2))  &
                           * (loc_x_node(2)+loc_x_node(3))  &
                           + (loc_z_node(1)-loc_z_node(3))  &
                           * (loc_x_node(3)+loc_x_node(1)) )
        grid % sz(s)= 0.5 * ((loc_x_node(2)-loc_x_node(1))  &
                           * (loc_y_node(2)+loc_y_node(1))  &
                           + (loc_x_node(3)-loc_x_node(2))  &
                           * (loc_y_node(2)+loc_y_node(3))  &
                           + (loc_x_node(1)-loc_x_node(3))  &
                           * (loc_y_node(3)+loc_y_node(1)) )
      else
        print *, '# Compute_Grid_Geometry: something horrible has happened !'
        stop
      end if

      ! Barycenters
      if(grid % faces_n_nodes(s) .eq. 4) then  
        grid % xf(s) = (   loc_x_node(1)+loc_x_node(2)          &
                         + loc_x_node(3)+loc_x_node(4) ) / 4.0
        grid % yf(s) = (   loc_y_node(1)+loc_y_node(2)          &
                         + loc_y_node(3)+loc_y_node(4) ) / 4.0
        grid % zf(s) = (   loc_z_node(1)+loc_z_node(2)          &
                         + loc_z_node(3)+loc_z_node(4) ) / 4.0
      else if(grid % faces_n_nodes(s) .eq. 3) then  
        grid % xf(s) = (loc_x_node(1)+loc_x_node(2)+loc_x_node(3)) / 3.0
        grid % yf(s) = (loc_y_node(1)+loc_y_node(2)+loc_y_node(3)) / 3.0
        grid % zf(s) = (loc_z_node(1)+loc_z_node(2)+loc_z_node(3)) / 3.0
      end if 

    end do ! through sides

    !--------------------------------------!
    !   Calculate boundary cell centers    !
    !--------------------------------------!
    !   => depends on: xc,yc,zc,Sx,Sy,Sz   !
    !   <= gives:      xc,yc,zc for c<0    !   
    !--------------------------------------!
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      tot_surf = sqrt(grid % sx(s)*grid % sx(s) +  &
                      grid % sy(s)*grid % sy(s) +  &
                      grid % sz(s)*grid % sz(s))

      if(c2  < 0) then
        t = (   grid % sx(s) * (grid % xf(s)-grid % xc(c1))           &
              + grid % sy(s) * (grid % yf(s)-grid % yc(c1))           &
              + grid % sz(s) * (grid % zf(s)-grid % zc(c1)) ) / tot_surf
        grid % xc(c2) = grid % xc(c1) + grid % sx(s)*t / tot_surf
        grid % yc(c2) = grid % yc(c1) + grid % sy(s)*t / tot_surf
        grid % zc(c2) = grid % zc(c1) + grid % sz(s)*t / tot_surf
      end if 
    end do ! through sides

    !---------------------------------------------!
    !   Find the sides on the periodic boundary   !
    !---------------------------------------------!
    !   => depends on: xc,yc,zc,Sx,Sy,Sz          !
    !   <= gives:      Dx,Dy,Dz                   !
    !---------------------------------------------!
    if(rrun) then
    n_per = 0
    do s = 1, grid % n_faces

      ! Initialize
      grid % dx(s) = 0.0
      grid % dy(s) = 0.0
      grid % dz(s) = 0.0

      c1 = grid % faces_c(1, s)
      c2 = grid % faces_c(2, s)
      if(c2   >  0) then

        ! Scalar product of the side with line c1-c2 is a good criterion
        if( (grid % sx(s) * (grid % xc(c2) - grid % xc(c1) ) +  &
             grid % sy(s) * (grid % yc(c2) - grid % yc(c1) ) +  &
             grid % sz(s) * (grid % zc(c2) - grid % zc(c1) ))  < 0.0 ) then

          n_per = n_per + 1
 
          ! Find the coordinates of ...
          m = face_c_to_c(s,2)

          if(grid % faces_n_nodes(s) .eq. 4) then   

            ! Coordinates of the shadow face
            xs2=.25*(  grid % xn(grid % cells_n(f4n(m,1), c2))  &
                     + grid % xn(grid % cells_n(f4n(m,2), c2))  &
                     + grid % xn(grid % cells_n(f4n(m,3), c2))  &
                     + grid % xn(grid % cells_n(f4n(m,4), c2)))

            ys2=.25*(  grid % yn(grid % cells_n(f4n(m,1), c2))  &
                     + grid % yn(grid % cells_n(f4n(m,2), c2))  &
                     + grid % yn(grid % cells_n(f4n(m,3), c2))  &
                     + grid % yn(grid % cells_n(f4n(m,4), c2)))

            zs2=.25*(  grid % zn(grid % cells_n(f4n(m,1), c2))  &
                     + grid % zn(grid % cells_n(f4n(m,2), c2))  &
                     + grid % zn(grid % cells_n(f4n(m,3), c2))  &
                     + grid % zn(grid % cells_n(f4n(m,4), c2)))
 
          else if(grid % faces_n_nodes(s) .eq. 3) then  

            ! Coordinates of the shadow face
            xs2 = ONE_THIRD * (grid % xn(grid % cells_n(f3n(m,1), c2))  &
                             + grid % xn(grid % cells_n(f3n(m,2), c2))  &
                             + grid % xn(grid % cells_n(f3n(m,3), c2)) )

            ys2 = ONE_THIRD * (grid % yn(grid % cells_n(f3n(m,1), c2))  &
                             + grid % yn(grid % cells_n(f3n(m,2), c2))  &
                             + grid % yn(grid % cells_n(f3n(m,3), c2)) )

            zs2 = ONE_THIRD * (grid % zn(grid % cells_n(f3n(m,1), c2))  &
                             + grid % zn(grid % cells_n(f3n(m,2), c2))  &
                             + grid % zn(grid % cells_n(f3n(m,3), c2)) )

          end if 

          grid % dx(s) = grid % xf(s) - xs2  !------------------------!
          grid % dy(s) = grid % yf(s) - ys2  ! later: xc2 = xc2 + Dx  !
          grid % dz(s) = grid % zf(s) - zs2  !------------------------!

        end if !  S*(c2-c1) < 0.0
      end if  !  c2 > 0
    end do    !  sides  
    print '(a38,i7)', '# Number of periodic faces:          ', n_per
    end if

  !----------------------------------!
  !   Calculate the cell volumes     !
  !----------------------------------!
  !   => depends on: xc,yc,zc,       !
  !                  Dx,Dy,Dz,       !
  !                  xsp, ysp, zsp   !
  !   <= gives:      volume          !
  !----------------------------------!
  if(rrun) then
    do c = 1, grid % n_cells
      grid % vol(c)=0.0
    end do

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)   

      do n = 1, grid % faces_n_nodes(s)  ! for quadrilateral an triangular faces
        loc_x_node(n) = grid % xn(grid % faces_n(n,s))
        loc_y_node(n) = grid % yn(grid % faces_n(n,s))
        loc_z_node(n) = grid % zn(grid % faces_n(n,s))
      end do   

      ! First cell
      x_cell_tmp = grid % xc(c1)
      y_cell_tmp = grid % yc(c1)
      z_cell_tmp = grid % zc(c1)
      dsc1 = Distance(x_cell_tmp,   y_cell_tmp,   z_cell_tmp,    &
                      grid % xf(s), grid % yf(s), grid % zf(s)) 
      grid % vol(c1) = grid % vol(c1)                                      &
                 + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),      &
                              loc_x_node(1),loc_y_node(1),loc_z_node(1),   &
                              loc_x_node(2),loc_y_node(2),loc_z_node(2),   &
                              x_cell_tmp,y_cell_tmp,z_cell_tmp)
      grid % vol(c1) = grid % vol(c1)                                      &
                 + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),      &
                              loc_x_node(2),loc_y_node(2),loc_z_node(2),   &
                              loc_x_node(3),loc_y_node(3),loc_z_node(3),   &
                              x_cell_tmp,y_cell_tmp,z_cell_tmp)
      if(grid % faces_n_nodes(s) .eq. 4) then
        grid % vol(c1) = grid % vol(c1)                                    &
                   + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),    &
                                loc_x_node(3),loc_y_node(3),loc_z_node(3), &
                                loc_x_node(4),loc_y_node(4),loc_z_node(4), &
                                x_cell_tmp,y_cell_tmp,z_cell_tmp)
        grid % vol(c1) = grid % vol(c1)                                    &
                   + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),    &
                                loc_x_node(4),loc_y_node(4),loc_z_node(4), &
                                loc_x_node(1),loc_y_node(1),loc_z_node(1), &
                                x_cell_tmp,y_cell_tmp,z_cell_tmp)
      else if(grid % faces_n_nodes(s) .eq. 3) then
        grid % vol(c1) = grid % vol(c1)                                    &
                   + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),    &
                                loc_x_node(3),loc_y_node(3),loc_z_node(3), &
                                loc_x_node(1),loc_y_node(1),loc_z_node(1), &
                                x_cell_tmp,y_cell_tmp,z_cell_tmp)
      end if

      ! Second cell
      if(c2  > 0) then
        x_cell_tmp = grid % xc(c2) + grid % dx(s)
        y_cell_tmp = grid % yc(c2) + grid % dy(s)
        z_cell_tmp = grid % zc(c2) + grid % dz(s)
        dsc2=Distance(x_cell_tmp,   y_cell_tmp,   z_cell_tmp,    &
                      grid % xf(s), grid % yf(s), grid % zf(s)) 
        grid % vol(c2) = grid % vol(c2)                                    &
                   - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),    &
                                loc_x_node(1),loc_y_node(1),loc_z_node(1), &
                                loc_x_node(2),loc_y_node(2),loc_z_node(2), &
                                x_cell_tmp,y_cell_tmp,z_cell_tmp)
        grid % vol(c2) = grid % vol(c2)                                    &
                   - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),    &
                                loc_x_node(2),loc_y_node(2),loc_z_node(2), &
                                loc_x_node(3),loc_y_node(3),loc_z_node(3), &
                                x_cell_tmp,y_cell_tmp,z_cell_tmp)
        if(grid % faces_n_nodes(s) .eq. 4) then
          grid % vol(c2) = grid % vol(c2)                                  &
                     - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),  &
                           loc_x_node(3),loc_y_node(3),loc_z_node(3),      &
                           loc_x_node(4),loc_y_node(4),loc_z_node(4),      &
                           x_cell_tmp,y_cell_tmp,z_cell_tmp)
          grid % vol(c2) = grid % vol(c2)                                  &
                     - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),  &
                           loc_x_node(4),loc_y_node(4),loc_z_node(4),      &
                           loc_x_node(1),loc_y_node(1),loc_z_node(1),      &
                           x_cell_tmp,y_cell_tmp,z_cell_tmp)
        else if(grid % faces_n_nodes(s) .eq. 3) then
          grid % vol(c2) = grid % vol(c2)                                  &
                     - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),  &
                           loc_x_node(3),loc_y_node(3),loc_z_node(3),      &
                           loc_x_node(1),loc_y_node(1),loc_z_node(1),      &
                           x_cell_tmp,y_cell_tmp,z_cell_tmp)
        end if
      else
        dsc2 = 0.0
      end if

    end do
  end if

  !------------------------------------------------------------!
  !   Calculate the interpolation factors for the cell sides   !
  !------------------------------------------------------------!
  if(rrun) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! First cell
      xc1 = grid % xc(c1)
      yc1 = grid % yc(c1)
      zc1 = grid % zc(c1)
      dsc1=Distance(xc1, yc1, zc1, grid % xf(s), grid % yf(s), grid % zf(s))

      ! Second cell (pls. check if xsi=xc on the boundary)
      xc2 = grid % xc(c2) + grid % dx(s)
      yc2 = grid % yc(c2) + grid % dy(s)
      zc2 = grid % zc(c2) + grid % dz(s)
      dsc2=Distance(xc2, yc2, zc2, grid % xf(s), grid % yf(s), grid % zf(s))

      ! Interpolation factor
      grid % f(s) = dsc2 / (dsc1 + dsc2)
    end do
  end if

  end subroutine
