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
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  type(Var_Type),      pointer :: u, v, w
  type(Particle_Type), pointer :: part
  logical,             pointer :: deposited            ! part. deposition flag
  logical,             pointer :: escaped              ! part. departure  flag
  integer                      :: c, c2                ! nearest cell
  real                         :: rx_nx_o, ry_ny_o, rz_nz_o,  &
                                  rx_nx_n, ry_ny_n, rz_nz_n,  &
                                  nx,     ny,     nz
!==============================================================================!

  ! Take aliases
  flow      => swarm % pnt_flow
  grid      => swarm % pnt_grid
  u         => flow % u
  v         => flow % v
  w         => flow % w
  part      => swarm % particle(k)
  deposited => part  % deposited
  escaped   => part  % escaped

  c  = part % cell      ! index of the closest cell for interpolation
  c2 = part % bnd_cell  ! index of the closest boundary cell for reflection

  ! Normal to the wall face
  nx = grid % sx(part % bnd_face) / grid % s(part % bnd_face)
  ny = grid % sy(part % bnd_face) / grid % s(part % bnd_face)
  nz = grid % sz(part % bnd_face) / grid % s(part % bnd_face)

  !-------------------------------------------------!
  !   Check if particle is approaching a boundary   !
  !-------------------------------------------------!
  if(   part % u * nx  &
      + part % v * ny  &
      + part % w * nz >= 0.0 ) then

    ! Old dot product of rx and nx
    rx_nx_o = (grid % xc(c2) - part % x_o) * nx
    ry_ny_o = (grid % yc(c2) - part % y_o) * ny
    rz_nz_o = (grid % zc(c2) - part % z_o) * nz

    ! New dot product of rx and nx
    rx_nx_n = (grid % xc(c2) - part % x) * nx
    ry_ny_n = (grid % yc(c2) - part % y) * ny
    rz_nz_n = (grid % zc(c2) - part % z) * nz

    if( (   rx_nx_o * rx_nx_n  &
          + ry_ny_o * ry_ny_n  &
          + rz_nz_o * rz_nz_n ) <= 0.0 ) then
      PRINT *, k, 'particle will hit the wall!'
      PRINT *, 'near wall face is: ', grid % cells_bnd_face(c2)
      PRINT *, 'boundary: ', grid % xc(c2), grid % yc(c2), grid % zc(c2)
      PRINT *, 'particle: ', part % x, part % y, part % z
      PRINT *, 'inside: ',   grid % xc(c), grid % yc(c), grid % zc(c)
      PRINT *, 'normal: ', nx, ny, nz
    end if

    !-----------------------------------!
    !    Trap condition (deposition)    !
    !-----------------------------------!
    if(swarm % rst <= TINY .and. .not. deposited) then
      if(part % y <= -0.005) then    !just for the moment
        deposited = .true.
        swarm % cnt_d = swarm % cnt_d + 1
        print *, k, 'Particle is deposited!'
      end if
    end if  ! trap condition

    !--------------------------!
    !   Reflection condition   !
    !--------------------------!
    if(swarm % rst > TINY) then
      if(part % y <= -0.005) then    !just for the moment until i make it generic
        part % y = -0.004999

        ! Change the direction of velocity
        part % u = part % u * ( swarm % rst)
        part % v = part % v * (-swarm % rst)
        part % w = part % w * ( swarm % rst)

        ! Update the particle position after reflection
        part % x = part % x + part % u * swarm % dt
        part % y = part % y + part % v * swarm % dt
        part % z = part % z + part % w * swarm % dt

        ! Increasing the number of particle reflections
        swarm % cnt_r = swarm % cnt_r + 1   ! to be engineered because ...
                                            ! ... a single particle can ...
                                            ! ... bounce several times.
        print *, k, 'Particle is reflected!'
      end if
    end if  ! reflection condition

    !--------------------------------------------------------------!
    !   Departure condition (particles escaping from the domain)   !
    !--------------------------------------------------------------!
    if(part % x .ge. 0.104999999999    &
      .or. part % x .le.  -0.104999999999) then    !just for the moment
      escaped =  .true.
      swarm % cnt_e = swarm % cnt_e + 1
      print *,k,'Particle escaped from outlet!'
    end if

  end if  ! approaching a boundary

  end subroutine
