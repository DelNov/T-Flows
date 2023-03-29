!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Front_Type), pointer :: Front
  type(Var_Type),   pointer :: fun, t, smooth
  type(Vof_Type),   pointer :: t_sat, latent_heat
  integer                   :: c, fu, s, e, nb, nc
  real                      :: t_s, x, y, z, r_x, r_y, r_z
  real                      :: radius, dist
  real                      :: dx, dy, dz, grad
  character(len=80)         :: filename  ! don't use SL for separate compilation
  integer                   :: n, cn, i_cel, i_nod, i_fac, reg
  real                      :: xe, ye, ze, xcn, ycn, zcn, dg
  real                      :: temp_b, grad_norm_x, grad_norm_y, grad_norm_z
  real                      :: xb, yb, zb, c_size
  real                      :: dotprod, dist_cn, dist_bcn
  real, dimension(3)        :: n_v, c_cent, b_v, cn_v
  integer, pointer, contiguous  :: cell_list(:)
  real, pointer, contiguous     :: cell_dist(:)
!==============================================================================!

  ! Take aliases
  Grid  => Flow % pnt_grid
  Front => Vof % Front
  fun   => Vof % fun
  smooth => Vof % smooth
  t     => Flow % t
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells
  
  t_s = Vof % t_sat
  r_x = 0.0
  r_y = 0.0
  r_z = 0.0
  radius = 0.002
  
!  do s = 1, Grid % n_bnd_cells
!    dx = abs(Grid % dx(s))
!    dy = abs(Grid % dy(s))
!    dz = abs(Grid % dz(s))
!  end do
!   Flow % u % n(:) =   0.0
!   Flow % v % n(:) =   0.0
  !---------------------------------!
  !   Initialize T values           !
  !---------------------------------!  

  ! For every cell within radius t=t_sat
  do reg = Boundary_And_Inside_Regions()
    do c = Cells_In_Region(reg)
      x = Grid % xc(c)
      y = Grid % yc(c)
      z = Grid % zc(c)
      ! Sphere
  !   dist = sqrt((x-r_x)**2 + (y-r_y)**2 + (z-r_z)**2)
      ! Cylinder
      dist = sqrt((x-r_x)**2 + (z-r_z)**2)

      if(fun % n(c) < 0.5) then
        if(Flow % heat_transfer) then
           Flow % t % n(c) = t_s 
           Flow % t % o(c) = t_s 
        end if 
      else 
        if(Flow % heat_transfer) then
           Flow % t % n(c) = t_s + 10000*(dist-radius)
           Flow % t % o(c) = t_s + 10000*(dist-radius)
        end if
      end if
    end do
  end do

  call Flow % Grad_Variable(t)
  ! Calculate smooth variable from vof function ...
  call Vof % Smooth_Scalar(Grid, fun % n,   &
                           smooth % n(-nb:nc), Vof % n_conv_curv)

  ! ... and find its gradients as well
  call Flow % Grad_Variable(smooth)

  
  if(DEBUG) then 
   call Work % Connect_Real_Cell(cell_dist)
   call Work % Connect_Int_Cell(cell_list)
   
   do c = 1, Grid % n_cells
      e = Vof % Front % elem_in_cell(c)
      if(e > 0) then
        x = Grid % xc(c)
        y = Grid % yc(c)
        z = Grid % zc(c)
        xe = Vof % Front % Elem(e) % xe
        ye = Vof % Front % Elem(e) % ye
        ze = Vof % Front % Elem(e) % ze

        ! Unity normal vector
        n_v = (/Vof % nx(c), Vof % ny(c), Vof % nz(c)/)

        ! Work out the cell size locally
        c_size = HUGE
        do i_fac = 1, Grid % cells_n_faces(c)
          dx = abs(Grid % dx(i_fac))
          dy = abs(Grid % dy(i_fac))
          dz = abs(Grid % dz(i_fac))
          c_size = min(c_size, sqrt((dx**2 + dy**2 + dz**2)))
        end do
        dg = 1.2 * c_size
        b_v = n_v * dg
        xb = xe + b_v(1)
        yb = ye + b_v(2)
        zb = ze + b_v(3)
        write(2000, *) xb, yb, zb
       ! Liquid phase    
        if(Vof % fun % n(c) > 0.5) then
          cell_dist(:) = HUGE
          cell_list(:) = 0

          ! Calculate the distances of all cells around
          ! the node and store them in cell_dist (and cell_list)
          do i_nod = 1, abs(Grid % cells_n_nodes(c)) ! local
            n = Grid % cells_n(i_nod, c)             ! global
            do i_cel = 1, Grid % nodes_n_cells(n)    ! local
              cn = Grid % nodes_c(i_cel, n)          ! global

              if(cn > 0                  .and.  &
                 cn .ne. c               .and.  &
                 Vof % fun % n(cn) > 0.5 .and.  &
                 Vof % Front % elem_in_cell(cn) .eq. 0) then
                xcn = Grid % xc(cn)
                ycn = Grid % yc(cn)
                zcn = Grid % zc(cn)

                cell_dist(i_cel) = sqrt( (xb - xcn)**2     &
                                       + (yb - ycn)**2     &
                                       + (zb - zcn)**2 )
                cell_list(i_cel) = cn
              end if
            end do
          end do

          ! Find out the cell with shortest distance to b
          call Sort % Real_Carry_Int(cell_dist(1:i_cel), cell_list(1:i_cel))

          ! Closest cell's index is now in cell_list(1)
          cn  = cell_list(1)
          xcn = Grid % xc(cn)
          ycn = Grid % yc(cn)
          zcn = Grid % zc(cn)
          write(3000, *) xcn, ycn, zcn
!          temp_b = Flow % t % n(cn) + (Flow % t % x(cn) * (xb - xcn) &
!                                     + Flow % t % y(cn) * (yb - ycn) &
!                                     + Flow % t % z(cn) * (zb - zcn))
                                     
          ! Calculate normal gradient in liquid phase
!          norm_grad_1(c) = (temp_b - Vof % t_sat) / dg
        end if
       
       ! Vapor phase
        if(Vof % fun % n(c) < 0.5) then
          b_v = -n_v * dg
          xb = xe + b_v(1)
          yb = ye + b_v(2)
          zb = ze + b_v(3)
          cell_dist(:) = HUGE
          cell_list(:) = 0
          write(4000, *) xb, yb, zb
          ! Calculate the distances of all cells around
          ! the node and store them in cell_dist (and cell_list)
          do i_nod = 1, abs(Grid % cells_n_nodes(c)) ! local
            n = Grid % cells_n(i_nod, c)             ! global
            do i_cel = 1, Grid % nodes_n_cells(n)    ! local
              cn = Grid % nodes_c(i_cel, n)          ! global

              if(cn .ne. c               .and.  &
                 Vof % fun % n(cn) < 0.5 .and.  &
                 Vof % Front % elem_in_cell(cn) .eq. 0) then
                xcn = Grid % xc(cn)
                ycn = Grid % yc(cn)
                zcn = Grid % zc(cn)

                cell_dist(i_cel) = sqrt(  (xb - xcn)**2  &
                                        + (yb - ycn)**2  &
                                        + (zb - zcn)**2 )
                cell_list(i_cel) = cn
              end if
            end do
          end do

          ! Find out the cell with shortest distance to b
          call Sort % Real_Carry_Int(cell_dist(1:i_cel), cell_list(1:i_cel))

          ! Closest cell's index is now in cell_list(1)
          cn  = cell_list(1)
          xcn = Grid % xc(cn)
          ycn = Grid % yc(cn)
          zcn = Grid % zc(cn)
          write(5000, *) xcn, ycn, zcn
!          temp_b = Flow % t % n(cn) + (Flow % t % x(cn) * (xb - xcn) &
!                                     + Flow % t % y(cn) * (yb - ycn) &
!                                     + Flow % t % z(cn) * (zb - zcn))
          ! Calculate normal gradient in vapor phase
!          norm_grad_0(c) = (temp_b - Vof % t_sat) / dg
        end if
      end if
    end do
    call Work % Disconnect_Real_Cell(cell_dist)
    call Work % Disconnect_Int_Cell(cell_list)
  end if
 
!          do i_nod = 1, abs(Grid % cells_n_nodes(c)) ! local
!            n = Grid % cells_n(i_nod, c)             ! global
!            do i_cel = 1, Grid % nodes_n_cells(n)    ! local
!              cn = Grid % nodes_c(i_cel, n)          ! global
!              if(cn > 0                                .and.  &
!                 cn .ne. c                             .and.  &
!                 Vof % Front % elem_in_cell(cn) .eq. 0 .and.  &
!                 Vof % fun % n(cn) > 0.5) then
!                xcn = Grid % xc(cn)
!                ycn = Grid % yc(cn)
!                zcn = Grid % zc(cn)
!                if((xcn - xp) > 0.0 .and. (ycn - yp) > 0.0 .and. &
!                   (zcn - zp) > 0.0) then
!                  db(cn) = sqrt((xcn - xp)**2 + (ycn - yp)**2 &
!                         + (zcn - zp)**2)
!                end if
!                min_d = minval(db)
!                if(min_d > 0.0) then
!                  PRINT*, "P", xp, yp, zp
!                  PRINT*, "CN", xcn, ycn, zcn
!                end if 
!                
!                d3 = (/xcn - xe, ycn - ye, zcn - ze/)
!                d3_n = dot_product(d3, d1_n/norm2(d1_n))*d1_n/norm2(d1_n)
!                stop
!                if(norm2(d3) > norm2(d1_n)) then
!                  write(filename,'("test_v2-",i7.7,".vtk")') cn  
!                  open(newunit=fu, file=filename)
!                  write(fu,'(a26)')     '# vtk DataFile Version 2.0'
!                  write(fu,'(a6,i7.7)') 'File: ', e
!                  write(fu,'(a5)')      'ASCII'
!                  write(fu,*)           ' '
!                  write(fu,'(a16)')     'DATASET POLYDATA'
!                  write(fu,'(a26)')     'POINTS   3   float'
!                  write(fu,'(3f12.6)') Front % Elem(e) % xe, &
!                                       Front % Elem(e) % ye, &
!                                       Front % Elem(e) % ze
!                  write(fu,'(3f12.6)') d1_n(1), &
!                                       d1_n(2), &
!                                       d1_n(3)
!                  write(fu,'(3f12.6)') d3_n(1), &
!                                       d3_n(2), &
!                                       d3_n(3)
!                  write(fu,'(a26)') 'POLYGONS     2     6'
!                  write(fu,'(a8)')  '    2'
!                  write(fu,'(a8)')  '    0'
!                  write(fu,'(a8)')  '    1'
!                  write(fu,'(a8)')  '    2'
!                  write(fu,'(a8)')  '    0'
!                  write(fu,'(a8)')  '    2' 
!                  close(fu)
!                end if
!                
!                
!              end if
!            end do
!          end do
!        end if
!      end if
!    end do

  end subroutine

