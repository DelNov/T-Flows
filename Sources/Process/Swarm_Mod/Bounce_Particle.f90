!==============================================================================!
  subroutine Swarm_Mod_Bounce_Particle(swarm, k)
!------------------------------------------------------------------------------!
!   Interpolates velocity at the particle's position                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer, intent(in)      :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  logical,             pointer :: deposited            ! part. deposition flag
  logical,             pointer :: escaped              ! part. departure  flag
  integer                      :: c, c2, s, bs         ! nearest cells, face
  real                         :: vel_dot_n
  real                         :: xi, yi, zi, nx, ny, nz
  real                         :: u_ref, v_ref, w_ref
  real                         :: part_x_n, part_y_n, part_z_n
  real                         :: part_x_o, part_y_o, part_z_o
  real                         :: vec_new_face_dot_n, vec_old_face_dot_n
  real                         :: lx, ly, lz, dsc_o, dsc_n
!==============================================================================!

  ! Take aliases
  grid      => swarm % pnt_grid
  part      => swarm % particle(k)
  deposited => part  % deposited
  escaped   => part  % escaped

  c  = part % cell      ! index of the closest cell for interpolation
  c2 = part % bnd_cell  ! index of the closest boundary cell for reflection
  s  = 0
  if(c2 .ne. 0) then
    s = grid % cells_bnd_face(c2)  ! index of the closest boundary face
  end if

  ! Debug: write(this_proc*1000+k, '(a,4i7)')  &
  ! Debug:                         'time, part, cell and face: ', n, k, c2, s

  ! If no boundary cell close, no need to stay here
  if(c2 .eq. 0) return

  !---------------------------------!
  !                                 !
  !                                 !
  !   Face at the wall is defined   !
  !                                 !
  !                                 !
  !---------------------------------!
  if(s .ne. 0) then

    ! Normal to the wall face
    ! (It points outwards, but you change the sign!!!)
    nx = -grid % sx(s) / grid % s(s)
    ny = -grid % sy(s) / grid % s(s)
    nz = -grid % sz(s) / grid % s(s)

    ! Velocity normal to the wall
    vel_dot_n = part % u * nx   &
              + part % v * ny   &
              + part % w * nz

    ! Vector connecting particle with boundary face center; new and old
    vec_new_face_dot_n = (part % x_n - grid % xf(s)) * nx  &
                       + (part % y_n - grid % yf(s)) * ny  &
                       + (part % z_n - grid % zf(s)) * nz

    vec_old_face_dot_n = (part % x_o - grid % xf(s)) * nx  &
                       + (part % y_o - grid % yf(s)) * ny  &
                       + (part % z_o - grid % zf(s)) * nz

    !---------------------------------------------------------------------!
    !                                                                     !
    !   Trouble!  Particle is out, but wasn't detected by the algorithm   !
    !   most probably because time step was too big, or because no sub-   !
    !   time stepping was used.                                           !
    !                                                                     !
    !---------------------------------------------------------------------!
    if( vec_new_face_dot_n < 0.0 .and. vec_old_face_dot_n < 0.0) then

      ! Try to move the particle back and see if it crosses
      do bs = 1, 16

        ! Move the particle back
        part_x_n = part % x_n - swarm % dt * part % u * bs
        part_y_n = part % y_n - swarm % dt * part % v * bs
        part_z_n = part % z_n - swarm % dt * part % w * bs
        part_x_o = part % x_o - swarm % dt * part % u * bs
        part_y_o = part % y_o - swarm % dt * part % v * bs
        part_z_o = part % z_o - swarm % dt * part % w * bs

        ! Vector connecting particle with boundary face center; new and old
        vec_new_face_dot_n = (part_x_n - grid % xf(s)) * nx  &
                           + (part_y_n - grid % yf(s)) * ny  &
                           + (part_z_n - grid % zf(s)) * nz

        vec_old_face_dot_n = (part_x_o - grid % xf(s)) * nx  &
                           + (part_y_o - grid % yf(s)) * ny  &
                           + (part_z_o - grid % zf(s)) * nz

        ! Managed to recover the particle
        if( vec_new_face_dot_n < 0.0 .and. vec_old_face_dot_n > 0.0) then
          print '(a)', ' #======================================' //  &
                       '=================='
          print '(a)', ' # WARNING!'
          print '(a,i7,a,i3,a)', ' # Particle ', k,  &
                                 ' got lost, but recovered in', bs, ' steps.'
          print '(a)', ' # Reducing the time step or increasing the number of'
          print '(a)', ' # sub time steps will hopefully help to resolve this'
          print '(a)', ' #-----------------------i--------------' //  &
                       '------------------'
          goto 1
        end if
      end do

      ! Failed to recovere the particle
      print '(a)', ' #======================================================='
      print '(a)', ' # PANIC!'
      print '(a,i7,a)', ' # Particle ', k, ' was lost and marked as "escaped"'
      print '(a)', ' # Reducing the time step or increasing the number of'
      print '(a)', ' # sub time steps will hopefully help to resolve this'
      print '(a)', ' #-------------------------------------------------------'
      escaped = .true.
      ! Mark the particle by enlarging it big time
      part % d = 1.0e-4
      ! Debug: write(this_proc*1000+k, '(a,2i7,a)')  &
      ! Debug:                         'time, part ', n, k, ' is out'
    end if

1   continue

    !--------------------------------------------------------------------!
    !                                                                    !
    !   Everything is fine, particle really crossed the boundary face    !
    !   or was recovered by backstepping (which isn't ideal of course)   !
    !                                                                    !
    !--------------------------------------------------------------------!
    if( vec_new_face_dot_n < 0.0 .and. vec_old_face_dot_n > 0.0) then

      ! Vector connecting old and new particle position
      lx = part % x_n - part % x_o
      ly = part % y_n - part % y_o
      lz = part % z_n - part % z_o

      ! Normalized distance from new to intersection
      dsc_n = (  (part % x_n - grid % xf(s)) * nx     &
               + (part % y_n - grid % yf(s)) * ny     &
               + (part % z_n - grid % zf(s)) * nz  )  &
            / (lx * nx + ly * ny + lz * nz)

      ! Normalized distance from old to intersection
      dsc_o = (  (grid % xf(s) - part % x_o) * nx     &
               + (grid % yf(s) - part % y_o) * ny     &
               + (grid % zf(s) - part % z_o) * nz  )  &
            / (lx * nx + ly * ny + lz * nz)

      ! Intersection point
      if(dsc_o > dsc_n) then
        xi = part % x_o + dsc_o * lx
        yi = part % y_o + dsc_o * ly
        zi = part % z_o + dsc_o * lz
      else
        xi = part % x_n - dsc_n * lx
        yi = part % y_n - dsc_n * ly
        zi = part % z_n - dsc_n * lz
      end if

      !---------------------------------!
      !   The boundary cell is a wall   !
      !---------------------------------!
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) == WALL .or.  & 
         Grid_Mod_Bnd_Cond_Type(grid, c2) == WALLFL) then

        ! Trap condition (deposition) >>> narrowed the tolerance  <<<
        if(swarm % rst <= TINY .or. abs(vel_dot_n) <= MILI) then

          deposited = .true.
          swarm % n_deposited(c2) = swarm % n_deposited(c2) + 1

          print '(a,i6,a,1pe11.3,1pe11.3,1pe11.3)',  &
                ' # Particle ', k, ' deposited at  : ', xi, yi, zi

          part % x_n = xi + part % d / 2.0 * nx
          part % y_n = yi + part % d / 2.0 * ny
          part % z_n = zi + part % d / 2.0 * nz

          part % x_o = part % x_n
          part % y_o = part % y_n
          part % z_o = part % z_n

          part % u = 0.0
          part % v = 0.0
          part % w = 0.0

        ! Reflected velocity
        else

          ! Work out reflected velocity
          u_ref = part % u - 2.0 * nx * vel_dot_n
          v_ref = part % v - 2.0 * ny * vel_dot_n
          w_ref = part % w - 2.0 * nz * vel_dot_n

          ! Set particle velocity to simply be the reflected (with restitution)
          part % u = u_ref * swarm % rst
          part % v = v_ref * swarm % rst
          part % w = w_ref * swarm % rst

          ! Place particle in reflected position
          part % x_n = xi + part % u * swarm % dt * dsc_n
          part % y_n = yi + part % v * swarm % dt * dsc_n
          part % z_n = zi + part % w * swarm % dt * dsc_n

          part % x_o = xi - part % u * swarm % dt * dsc_o
          part % y_o = yi - part % v * swarm % dt * dsc_o
          part % z_o = zi - part % w * swarm % dt * dsc_o

          ! Increasing the number of particle reflections
          swarm % n_reflected(c2) = swarm % n_reflected(c2) + 1

          ! Debug: write(this_proc*1000+k, '(a,2i7, 5(e12.4))')       &
          ! Debug:                         'time, part bounced at ',  &
          ! Debug:                         n, k,                      &
          ! Debug:                         xi, yi, zi, dsc_n, dsc_o

        end if  ! reflected velocity

      end if  ! boundary cell is a wall cell

      !------------------------------------!
      !   The boundary cell is an outlet   !
      !------------------------------------!
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) == OUTFLOW .or.  &
         Grid_Mod_Bnd_Cond_Type(grid, c2) == CONVECT) then
        escaped = .true.
        swarm % n_escaped(c2) = swarm % n_escaped(c2) + 1
      end if  ! it is an outflow

    end if  ! crossed a boundary cell

  end if  ! s .ne. 0

  end subroutine
