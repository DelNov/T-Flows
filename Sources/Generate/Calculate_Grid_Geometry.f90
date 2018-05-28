!==============================================================================!
  subroutine Calculate_Grid_Geometry(grid, rrun)
!------------------------------------------------------------------------------!
!   Calculates geometrical quantities of the grid.                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use gen_mod
  use Grid_Mod
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)     :: grid
  logical, intent(in) :: rrun     
!----------------------------------[Calling]-----------------------------------!
  real :: Tet_Volume        
  real :: Distance       
  real :: Distance_Squared       
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, m, n, s, new_face_1, new_face_2
  integer              :: b, n_wall_colors
  real                 :: local_x_node(4), local_y_node(4), local_z_node(4)
  real                 :: x_cell_tmp, y_cell_tmp, z_cell_tmp    
  real                 :: xs2, ys2, zs2
  real                 :: dsc1, dsc2          !  for the interpolation factors
  real                 :: t, tot_surf 
  real                 :: xc1, yc1, zc1, xc2, yc2, zc2 
  real                 :: x_min, x_max, y_min, y_max, z_min, z_max
  integer              :: f4n(6,4)
  integer              :: f3n(4,3)
  integer, allocatable :: wall_colors(:)
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
        local_x_node(n) = grid % xn(grid % faces_n(n,s))
        local_y_node(n) = grid % yn(grid % faces_n(n,s))
        local_z_node(n) = grid % zn(grid % faces_n(n,s))
      end do                       

      ! Cell side components
      if( grid % faces_n_nodes(s) .eq. 4 ) then
        grid % sx(s)= 0.5 * ((local_y_node(2)-local_y_node(1))  &
                           * (local_z_node(2)+local_z_node(1))  &
                           + (local_y_node(3)-local_y_node(2))  &
                           * (local_z_node(2)+local_z_node(3))  &
                           + (local_y_node(4)-local_y_node(3))  &
                           * (local_z_node(3)+local_z_node(4))  &
                           + (local_y_node(1)-local_y_node(4))  &
                           * (local_z_node(4)+local_z_node(1)) )
        grid % sy(s)= 0.5 * ((local_z_node(2)-local_z_node(1))  &
                           * (local_x_node(2)+local_x_node(1))  &
                           + (local_z_node(3)-local_z_node(2))  &
                           * (local_x_node(2)+local_x_node(3))  &
                           + (local_z_node(4)-local_z_node(3))  &
                           * (local_x_node(3)+local_x_node(4))  &
                           + (local_z_node(1)-local_z_node(4))  &
                           * (local_x_node(4)+local_x_node(1)) )
        grid % sz(s)= 0.5 * ((local_x_node(2)-local_x_node(1))  & 
                           * (local_y_node(2)+local_y_node(1))  & 
                           + (local_x_node(3)-local_x_node(2))  & 
                           * (local_y_node(2)+local_y_node(3))  &
                           + (local_x_node(4)-local_x_node(3))  & 
                           * (local_y_node(3)+local_y_node(4))  &
                           + (local_x_node(1)-local_x_node(4))  & 
                           * (local_y_node(4)+local_y_node(1)) )
      else if( grid % faces_n_nodes(s) .eq. 3 ) then 
        grid % sx(s)= 0.5 * ((local_y_node(2)-local_y_node(1))  &
                           * (local_z_node(2)+local_z_node(1))  & 
                           + (local_y_node(3)-local_y_node(2))  &
                           * (local_z_node(2)+local_z_node(3))  &
                           + (local_y_node(1)-local_y_node(3))  &
                           * (local_z_node(3)+local_z_node(1)) )
        grid % sy(s)= 0.5 * ((local_z_node(2)-local_z_node(1))  &
                           * (local_x_node(2)+local_x_node(1))  &
                           + (local_z_node(3)-local_z_node(2))  &
                           * (local_x_node(2)+local_x_node(3))  & 
                           + (local_z_node(1)-local_z_node(3))  &
                           * (local_x_node(3)+local_x_node(1)) )
        grid % sz(s)= 0.5 * ((local_x_node(2)-local_x_node(1))  &
                           * (local_y_node(2)+local_y_node(1))  &
                           + (local_x_node(3)-local_x_node(2))  & 
                           * (local_y_node(2)+local_y_node(3))  & 
                           + (local_x_node(1)-local_x_node(3))  & 
                           * (local_y_node(3)+local_y_node(1)) )
      else
        print *, '# Compute_Grid_Geometry: something horrible has happened !'
        stop
      end if

      ! Barycenters
      if(grid % faces_n_nodes(s) .eq. 4) then  
        grid % xf(s) = (   local_x_node(1)+local_x_node(2)          &
                         + local_x_node(3)+local_x_node(4) ) / 4.0
        grid % yf(s) = (   local_y_node(1)+local_y_node(2)          &
                         + local_y_node(3)+local_y_node(4) ) / 4.0
        grid % zf(s) = (   local_z_node(1)+local_z_node(2)          &
                         + local_z_node(3)+local_z_node(4) ) / 4.0
      else if(grid % faces_n_nodes(s) .eq. 3) then  
        grid % xf(s) = (local_x_node(1)+local_x_node(2)+local_x_node(3)) / 3.0
        grid % yf(s) = (local_y_node(1)+local_y_node(2)+local_y_node(3)) / 3.0
        grid % zf(s) = (local_z_node(1)+local_z_node(2)+local_z_node(3)) / 3.0
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
      endif 
    end do ! through sides

    !---------------------------------------------!
    !   Find the sides on the periodic boundary   !
    !---------------------------------------------!
    !   => depends on: xc,yc,zc,Sx,Sy,Sz          !
    !   <= gives:      Dx,Dy,Dz                   !
    !---------------------------------------------!
    if(rrun) then
    grid % n_sh = 0
    do s = 1, grid % n_faces

      ! Initialize
      grid % dx(s)=0.0
      grid % dy(s)=0.0
      grid % dz(s)=0.0

      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2   >  0) then

        ! Scalar product of the side with line c1-c2 is a good criterion
        if( (grid % sx(s) * (grid % xc(c2) - grid % xc(c1) ) +  &
             grid % sy(s) * (grid % yc(c2) - grid % yc(c1) ) +  &
             grid % sz(s) * (grid % zc(c2) - grid % zc(c1) ))  < 0.0 ) then

          grid % n_sh = grid % n_sh + 2
 
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
 
            ! Add shadow faces
            new_face_1 = grid % n_faces+grid % n_sh-1
            new_face_2 = grid % n_faces+grid % n_sh  
            grid % faces_n_nodes(new_face_1) = 4
            grid % faces_c(1, new_face_1) = c1 
            grid % faces_c(2, new_face_1) = -grid % n_bnd_cells-1
            grid % faces_n(1, new_face_1) = grid % faces_n(1,s)
            grid % faces_n(2, new_face_1) = grid % faces_n(2,s)
            grid % faces_n(3, new_face_1) = grid % faces_n(3,s)
            grid % faces_n(4, new_face_1) = grid % faces_n(4,s)
            grid % sx(new_face_1) = grid % sx(s)
            grid % sy(new_face_1) = grid % sy(s)
            grid % sz(new_face_1) = grid % sz(s)
            grid % xf(new_face_1) = grid % xf(s)
            grid % yf(new_face_1) = grid % yf(s)
            grid % zf(new_face_1) = grid % zf(s)
            grid % faces_n_nodes(new_face_2) = 4
            grid % faces_c(1, new_face_2) = c2 
            grid % faces_c(2, new_face_2) = -grid % n_bnd_cells-1
            grid % faces_n(1, new_face_2) = grid % cells_n(f4n(m,1), c2)
            grid % faces_n(2, new_face_2) = grid % cells_n(f4n(m,2), c2)
            grid % faces_n(3, new_face_2) = grid % cells_n(f4n(m,3), c2)
            grid % faces_n(4, new_face_2) = grid % cells_n(f4n(m,4), c2)
            grid % sx(new_face_2) = grid % sx(s)
            grid % sy(new_face_2) = grid % sy(s)
            grid % sz(new_face_2) = grid % sz(s)
            grid % xf(new_face_2) = xs2
            grid % yf(new_face_2) = ys2
            grid % zf(new_face_2) = zs2
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

            ! Add shadow faces
            new_face_1 = grid % n_faces+grid % n_sh-1
            new_face_2 = grid % n_faces+grid % n_sh  
            grid % faces_n_nodes(new_face_1) = 3
            grid % faces_c(1, new_face_1) = c1 
            grid % faces_c(2, new_face_1) = -grid % n_bnd_cells-1
            grid % faces_n(1, new_face_1) = grid % faces_n(1,s)
            grid % faces_n(2, new_face_1) = grid % faces_n(2,s)
            grid % faces_n(3, new_face_1) = grid % faces_n(3,s)
            grid % sx(new_face_1) = grid % sx(s)
            grid % sy(new_face_1) = grid % sy(s)
            grid % sz(new_face_1) = grid % sz(s)
            grid % xf(new_face_1) = grid % xf(s)
            grid % yf(new_face_1) = grid % yf(s)
            grid % zf(new_face_1) = grid % zf(s)
            grid % faces_n_nodes(new_face_2) = 3
            grid % faces_c(1, new_face_2) = c2 
            grid % faces_c(2, new_face_2) = -grid % n_bnd_cells-1
            grid % faces_n(1, new_face_2) = grid % cells_n(f3n(m,1), c2)
            grid % faces_n(2, new_face_2) = grid % cells_n(f3n(m,2), c2)
            grid % faces_n(3, new_face_2) = grid % cells_n(f3n(m,3), c2)
            grid % sx(new_face_2) = grid % sx(s)
            grid % sy(new_face_2) = grid % sy(s)
            grid % sz(new_face_2) = grid % sz(s)
            grid % xf(new_face_2) = xs2
            grid % yf(new_face_2) = ys2
            grid % zf(new_face_2) = zs2
          end if 

          grid % dx(s)=grid % xf(s)-xs2  !------------------------!
          grid % dy(s)=grid % yf(s)-ys2  ! later: xc2 = xc2 + Dx  !
          grid % dz(s)=grid % zf(s)-zs2  !------------------------!

        endif !  S*(c2-c1) < 0.0
      end if  !  c2 > 0
    end do    !  sides  
    print '(a38,i7)', '# Number of shadow faces:            ', grid % n_sh
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
      local_x_node(n) = grid % xn(grid % faces_n(n,s))
      local_y_node(n) = grid % yn(grid % faces_n(n,s))
      local_z_node(n) = grid % zn(grid % faces_n(n,s))
    end do   

    ! First cell
    x_cell_tmp = grid % xc(c1)
    y_cell_tmp = grid % yc(c1)
    z_cell_tmp = grid % zc(c1)
    dsc1 = Distance(x_cell_tmp,   y_cell_tmp,   z_cell_tmp,    &
                    grid % xf(s), grid % yf(s), grid % zf(s)) 
    grid % vol(c1) = grid % vol(c1)                                            &
               + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),            &
                            local_x_node(1),local_y_node(1),local_z_node(1),   &
                            local_x_node(2),local_y_node(2),local_z_node(2),   &
                            x_cell_tmp,y_cell_tmp,z_cell_tmp)
    grid % vol(c1) = grid % vol(c1)                                            &
               + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),            &
                            local_x_node(2),local_y_node(2),local_z_node(2),   &
                            local_x_node(3),local_y_node(3),local_z_node(3),   &
                            x_cell_tmp,y_cell_tmp,z_cell_tmp)
    if(grid % faces_n_nodes(s) .eq. 4) then
      grid % vol(c1) = grid % vol(c1)                                          &
                 + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),          &
                              local_x_node(3),local_y_node(3),local_z_node(3), &
                              local_x_node(4),local_y_node(4),local_z_node(4), &
                              x_cell_tmp,y_cell_tmp,z_cell_tmp)
      grid % vol(c1) = grid % vol(c1)                                          &
                 + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),          &
                              local_x_node(4),local_y_node(4),local_z_node(4), &
                              local_x_node(1),local_y_node(1),local_z_node(1), &
                              x_cell_tmp,y_cell_tmp,z_cell_tmp)
    else if(grid % faces_n_nodes(s) .eq. 3) then
      grid % vol(c1) = grid % vol(c1)                                          &
                 + Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),          &
                              local_x_node(3),local_y_node(3),local_z_node(3), &
                              local_x_node(1),local_y_node(1),local_z_node(1), &
                              x_cell_tmp,y_cell_tmp,z_cell_tmp)
    end if

    ! Second cell
    if(c2  > 0) then
      x_cell_tmp = grid % xc(c2) + grid % dx(s)
      y_cell_tmp = grid % yc(c2) + grid % dy(s)
      z_cell_tmp = grid % zc(c2) + grid % dz(s)
      dsc2=Distance(x_cell_tmp,   y_cell_tmp,   z_cell_tmp,    &
                    grid % xf(s), grid % yf(s), grid % zf(s)) 
      grid % vol(c2) = grid % vol(c2)                                          &
                 - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),          &
                              local_x_node(1),local_y_node(1),local_z_node(1), &
                              local_x_node(2),local_y_node(2),local_z_node(2), &
                              x_cell_tmp,y_cell_tmp,z_cell_tmp)
      grid % vol(c2) = grid % vol(c2)                                          &
                 - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),          &
                              local_x_node(2),local_y_node(2),local_z_node(2), &
                              local_x_node(3),local_y_node(3),local_z_node(3), &
                              x_cell_tmp,y_cell_tmp,z_cell_tmp)
      if(grid % faces_n_nodes(s) .eq. 4) then
        grid % vol(c2) = grid % vol(c2)                                    &
                   - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),    &
                         local_x_node(3),local_y_node(3),local_z_node(3),  &
                         local_x_node(4),local_y_node(4),local_z_node(4),  &
                         x_cell_tmp,y_cell_tmp,z_cell_tmp)
        grid % vol(c2) = grid % vol(c2)                                    &
                   - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),    &
                         local_x_node(4),local_y_node(4),local_z_node(4),  &
                         local_x_node(1),local_y_node(1),local_z_node(1),  &
                         x_cell_tmp,y_cell_tmp,z_cell_tmp)
      else if(grid % faces_n_nodes(s) .eq. 3) then
        grid % vol(c2) = grid % vol(c2)                                    &
                   - Tet_Volume(grid % xf(s),grid % yf(s),grid % zf(s),    &
                         local_x_node(3),local_y_node(3),local_z_node(3),  &
                         local_x_node(1),local_y_node(1),local_z_node(1),  &
                         x_cell_tmp,y_cell_tmp,z_cell_tmp)
      end if  
    else        
      dsc2=0.0
    end if

  end do
  end if

  !----------------------------------------------------------!
  !   Calculate:                                             ! 
  !      distance from the cell center to the nearest wall   !
  !----------------------------------------------------------!
  !   => depends on: xc,yc,zc inside and on the boundary     !
  !   <= gives:      wall_dist                               !
  !----------------------------------------------------------!
  grid % wall_dist = HUGE 

  call Grid_Mod_Print_Bnd_Cond_List(grid)
  print *, '#================================================================'
  if(rrun) then
    print *, '# Computing the distance to the walls (2/2)'           
  else            
    print *, '# Computing the distance to the walls (1/2)'           
  end if 
  print *, '# Type the list of boundary colors which represent walls,        '
  print *, '# separated by spaces.  These will be used for computation       '
  print *, '# of distance to the wall needed by some turbulence models.      '
  print *, '#----------------------------------------------------------------'
  call Tokenizer_Mod_Read_Line(5)
  n_wall_colors = line % n_tokens
  allocate(wall_colors(n_wall_colors))
  do b = 1, n_wall_colors
    read(line % tokens(b), *) wall_colors(b)
  end do
 
  if( (n_wall_colors .eq. 1) .and. (wall_colors(1) .eq. 0) ) then
    grid % wall_dist = 1.0
    print *, '# Distance to the wall set to 1.0 everywhere !'
  else 
    do b = 1, n_wall_colors
      do c1=1, grid % n_cells 
        do s = WallFacFst, WallFacLst      ! 1, grid % n_faces
          c2 = grid % faces_c(2,s)
          if(c2 < 0) then
            if(grid % bnd_cond % color(c2) .eq. wall_colors(b)) then
              grid % wall_dist(c1)=min(grid % wall_dist(c1), &
              Distance_Squared(grid % xc(c1),  &
                               grid % yc(c1),  &
                               grid % zc(c1),  &
                               grid % xf(s),   &
                               grid % yf(s),   &
                               grid % zf(s)))
            end if
          end if 
        end do
      end do
    end do

    do c = 1, grid % n_cells
      grid % wall_dist(c) = sqrt(grid % wall_dist(c))
    end do

    print *, '# Maximal distance to the wall: ',  &
                  maxval(grid % wall_dist(1:grid % n_cells))
    print *, '# Minimal distance to the wall: ',  &
                  minval(grid % wall_dist(1:grid % n_cells))
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
