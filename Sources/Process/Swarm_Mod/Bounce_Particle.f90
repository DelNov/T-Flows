!==============================================================================!
  subroutine Swarm_Mod_Bounce_Particle(swarm, k)
!------------------------------------------------------------------------------!
!   Interpolates velocity at the particle's position                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  logical,             pointer :: deposited            ! part. deposition flag
  logical,             pointer :: escaped              ! part. departure  flag
  integer                      :: c, c2, s             ! nearest cells, face
  real                         :: rx_nx_o, ry_ny_o, rz_nz_o,  &
                                  rx_nx_n, ry_ny_n, rz_nz_n,  &
                                  xi, yi, zi, dx, dy, dz, nx, ny, nz, f
!==============================================================================!

  ! Take aliases
  grid      => swarm % pnt_grid
  part      => swarm % particle(k)
  deposited => part  % deposited
  escaped   => part  % escaped

  c  = part % cell      ! index of the closest cell for interpolation
  c2 = part % bnd_cell  ! index of the closest boundary cell for reflection
  s  = part % bnd_face  ! index of the closest boundary face

  ! Normal to the wall face (it points outwards
  nx = grid % sx(s) / grid % s(s)
  ny = grid % sy(s) / grid % s(s)
  nz = grid % sz(s) / grid % s(s)

  !-------------------------------------------!
  !                                           !
  !                                           !
  !   If particle is approaching a boundary   !
  !                                           !
  !                                           !
  !-------------------------------------------!
  if(   part % u * nx  &
      + part % v * ny  &
      + part % w * nz >= 0.0 ) then

    ! Old dot product of rx and nx
    rx_nx_o = (grid % xc(c2) - part % x_o) * nx
    ry_ny_o = (grid % yc(c2) - part % y_o) * ny
    rz_nz_o = (grid % zc(c2) - part % z_o) * nz

    ! New dot product of rx and nx
    rx_nx_n = (grid % xc(c2) - part % x_n) * nx
    ry_ny_n = (grid % yc(c2) - part % y_n) * ny
    rz_nz_n = (grid % zc(c2) - part % z_n) * nz

    !------------------------------------------!
    !                                          !
    !   Did particle pass through a boundary   !
    !                                          !
    !------------------------------------------!
    if( (   rx_nx_o * rx_nx_n  &
          + ry_ny_o * ry_ny_n  &
          + rz_nz_o * rz_nz_n ) <= 0.0 ) then

      !--------------------------------------!
      !   Find where particle hit the wall   !
      !   (call this point i, like impact)   !
      !--------------------------------------!

      !---------------------------------------------!
      !   Equation for line which connects          !
      !   old and new particle position reads:      !
      !      ->        ->   ->                      !
      !   L: po + f * (pn - po)                     !
      !                                             !
      !   or:                                       !
      !   x(f) = xo + f * dx;  where dx = xn - xo   !
      !   y(f) = yo + f * dy;  where dy = yn - yo   !
      !   z(f) = zo + f * dz;  where dz = zn - zo   !
      !---------------------------------------------!
      dx = part % x_n - part % x_o
      dy = part % y_n - part % y_o
      dz = part % z_n - part % z_o

      !----------------------------------------------------!
      !   Vector connecting crossing and bondary           !
      !   center are orthogonal to face normal:            !
      !    ->   ->    ->                                   !
      !   (xi - xb) * n  = 0                               !
      !                                                    !
      !   Furhtermore, xi is on the line                   !
      !   connecting po and pn, hence:                     !
      !   ->    ->        ->                               !
      !   xi = po_x + f * dx                               !
      !                                                    !
      !   So the condition above reads                     !
      !     (po_x +  f * dx - xb) * nx                     !
      !   + (po_y +  f * dy - yb) * ny                     !
      !   + (po_z +  f * dz - zb) * nz = 0                 !
      !                                                    !
      !   And f is, eventually:                            !
      !                                                    !
      !       (xb-po_x)*nx + (yb-po_y)*ny + (zb-po_z)*nz   !
      !   f = ------------------------------------------   !
      !               dx*nx + dy*ny + dz*nz                !
      !                                                    !
      !----------------------------------------------------!
      f = (rx_nx_o + ry_ny_o + rz_nz_o) / (dx*nx + dy*ny + dz*nz)

      xi = part % x_o + f*dx
      yi = part % y_o + f*dy
      zi = part % z_o + f*dz

      !---------------------------------!
      !   The boundary cell is a wall   !
      !---------------------------------!
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) == WALL) then

        ! Trap condition (deposition)
        if(swarm % rst <= TINY) then
          deposited = .true.
          swarm % cnt_d = swarm % cnt_d + 1
          print *, k, 'Particle is deposited at: ', xi, yi, zi
        end if  ! trap condition

        ! Reflection condition   !
        if(swarm % rst > TINY) then
          part % y_n = -0.004999

          ! Change the direction of velocity
          part % u = part % u * ( swarm % rst)
          part % v = part % v * (-swarm % rst)
          part % w = part % w * ( swarm % rst)

          ! Update the particle position after reflection
          part % x_n = part % x_n + part % u * swarm % dt
          part % y_n = part % y_n + part % v * swarm % dt
          part % z_n = part % z_n + part % w * swarm % dt

          ! Increasing the number of particle reflections
          swarm % cnt_r = swarm % cnt_r + 1   ! to be engineered because ...
                                              ! ... a single particle can ...
                                              ! ... bounce several times.
          print *, k, 'Particle is reflected from: ', xi, yi, zi
        end if  ! reflection condition

      end if  ! is it a wall

      !------------------------------------!
      !   The boundary cell is an outlet   !
      !------------------------------------!
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) == OUTFLOW) then
        escaped =  .true.
        swarm % cnt_e = swarm % cnt_e + 1
        print *, k, 'Particle escaped from outlet at: ', xi, yi, zi
      end if  ! it is an outflow

    end if  ! really crossed a boundary cell

  end if  ! approaching a boundary

  end subroutine
