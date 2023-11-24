!==============================================================================!
  subroutine Bounce_Particle(Swarm, k)
!------------------------------------------------------------------------------!
!> This subroutine manages the interaction of particles with the boundaries of
!> the computational domain in T-Flows. It handles the reflection, deposition,
!> and escape of particles upon encountering walls, outflows, and other
!> boundary types.
!------------------------------------------------------------------------------!
! Functionality                                                                !
!                                                                              !
! * Particle-boundary Interaction: Detects and processes collisions between    !
!   particles and boundary surfaces.                                           !
! * Reflection and deposition: Determines whether a particle is reflected or   !
!   deposited based on its velocity and the nature of the boundary surface.    !
! * Escape handling: Identifies particles that exit the computational domain   !
!   through outflow boundaries and marks them as escaped.                      !
! * Position adjustment: Updates particle positions post-interaction to        !
!   reflect their new state (e.g., reflected position or deposited location).  !
! * Velocity update: Modifies particle velocities after reflection, taking     !
!   into account factors like the coefficient of restitution.                  !
! * Boundary crossing Validation: Ensures accurate detection of particles      !
!   crossing boundary faces, employing a backstepping method if necessary.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm of particles
  integer, intent(in)       :: k      !! particle number (rank)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Particle_Type), pointer :: Part
  logical,             pointer :: deposited            ! part. deposition flag
  logical,             pointer :: escaped              ! part. departure  flag
  integer                      :: c, c2, s, bs         ! nearest cells, face
  real                         :: vel_dot_n, vel_magn
  real                         :: xi, yi, zi, nx, ny, nz
  real                         :: u_ref, v_ref, w_ref
  real                         :: part_x_n, part_y_n, part_z_n
  real                         :: part_x_o, part_y_o, part_z_o
  real                         :: vec_new_face_dot_n, vec_old_face_dot_n
  real                         :: lx, ly, lz, dsc_o, dsc_n
!==============================================================================!

  ! Take aliases
  Grid      => Swarm % pnt_grid
  Part      => Swarm % Particle(k)
  deposited => Part  % deposited
  escaped   => Part  % escaped

  c  = Part % cell      ! index of the closest cell for interpolation
  c2 = Part % bnd_cell  ! index of the closest boundary cell for reflection
  s  = 0
  if(c2 .ne. 0) then
    s = Grid % cells_bnd_face(c2)  ! index of the closest boundary face
  end if

  ! Debug: write(This_Proc()*1000+k, '(a,4i7)')  &
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
    nx = -Grid % sx(s) / Grid % s(s)
    ny = -Grid % sy(s) / Grid % s(s)
    nz = -Grid % sz(s) / Grid % s(s)

    ! Vector connecting particle with boundary face center; new and old
    vec_new_face_dot_n = (Part % x_n - Grid % xf(s)) * nx  &
                       + (Part % y_n - Grid % yf(s)) * ny  &
                       + (Part % z_n - Grid % zf(s)) * nz

    vec_old_face_dot_n = (Part % x_o - Grid % xf(s)) * nx  &
                       + (Part % y_o - Grid % yf(s)) * ny  &
                       + (Part % z_o - Grid % zf(s)) * nz

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
        part_x_n = Part % x_n - Swarm % dt * Part % u * real(bs)
        part_y_n = Part % y_n - Swarm % dt * Part % v * real(bs)
        part_z_n = Part % z_n - Swarm % dt * Part % w * real(bs)
        part_x_o = Part % x_o - Swarm % dt * Part % u * real(bs)
        part_y_o = Part % y_o - Swarm % dt * Part % v * real(bs)
        part_z_o = Part % z_o - Swarm % dt * Part % w * real(bs)

        ! Vector connecting particle with boundary face center; new and old
        vec_new_face_dot_n = (part_x_n - Grid % xf(s)) * nx  &
                           + (part_y_n - Grid % yf(s)) * ny  &
                           + (part_z_n - Grid % zf(s)) * nz

        vec_old_face_dot_n = (part_x_o - Grid % xf(s)) * nx  &
                           + (part_y_o - Grid % yf(s)) * ny  &
                           + (part_z_o - Grid % zf(s)) * nz

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
      Part % d = 1.0e-4
      ! Debug: write(This_Proc()*1000+k, '(a,2i7,a)')  &
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
      lx = Part % x_n - Part % x_o
      ly = Part % y_n - Part % y_o
      lz = Part % z_n - Part % z_o

      ! Normalized distance from new to intersection
      dsc_n = (  (Part % x_n - Grid % xf(s)) * nx     &
               + (Part % y_n - Grid % yf(s)) * ny     &
               + (Part % z_n - Grid % zf(s)) * nz  )  &
            / (lx * nx + ly * ny + lz * nz)

      ! Normalized distance from old to intersection
      dsc_o = (  (Grid % xf(s) - Part % x_o) * nx     &
               + (Grid % yf(s) - Part % y_o) * ny     &
               + (Grid % zf(s) - Part % z_o) * nz  )  &
            / (lx * nx + ly * ny + lz * nz)

      ! Intersection point
      if(dsc_o > dsc_n) then
        xi = Part % x_o + dsc_o * lx
        yi = Part % y_o + dsc_o * ly
        zi = Part % z_o + dsc_o * lz
      else
        xi = Part % x_n - dsc_n * lx
        yi = Part % y_n - dsc_n * ly
        zi = Part % z_n - dsc_n * lz
      end if

      !---------------------------------!
      !   The boundary cell is a wall   !
      !---------------------------------!
      if(Grid % Bnd_Cond_Type(c2) == WALL .or.  & 
         Grid % Bnd_Cond_Type(c2) == WALLFL) then

        ! Velocity normal to the wall and velocity magnitude
        vel_dot_n = Part % u * nx   &
                  + Part % v * ny   &
                  + Part % w * nz
        vel_magn = sqrt(  Part % u * Part % u   &
                        + Part % v * Part % v   &
                        + Part % w * Part % w)

        ! Trap condition (deposition)
        if(Swarm % rst <= TINY                            .or.  &  ! sticky
           Swarm % rst <= 1.0-TINY .and. vel_magn < MICRO) then    ! too slow

          deposited = .true.
          Swarm % n_deposited(c2) = Swarm % n_deposited(c2) + 1

          print '(a,i6,a,1pe11.3,1pe11.3,1pe11.3)',  &
                ' # Particle ', k, ' deposited at  : ', xi, yi, zi

          Part % x_n = xi + Part % d / 2.0 * nx
          Part % y_n = yi + Part % d / 2.0 * ny
          Part % z_n = zi + Part % d / 2.0 * nz

          Part % x_o = Part % x_n
          Part % y_o = Part % y_n
          Part % z_o = Part % z_n

          Part % u = 0.0
          Part % v = 0.0
          Part % w = 0.0

        ! Reflected velocity
        else

          ! Work out reflected velocity
          u_ref = Part % u - 2.0 * nx * vel_dot_n
          v_ref = Part % v - 2.0 * ny * vel_dot_n
          w_ref = Part % w - 2.0 * nz * vel_dot_n

          ! Set particle velocity to simply be the reflected (with restitution)
          Part % u = u_ref * Swarm % rst
          Part % v = v_ref * Swarm % rst
          Part % w = w_ref * Swarm % rst

          ! Place particle in reflected position
          Part % x_n = xi + Part % u * Swarm % dt * dsc_n
          Part % y_n = yi + Part % v * Swarm % dt * dsc_n
          Part % z_n = zi + Part % w * Swarm % dt * dsc_n

          Part % x_o = xi - Part % u * Swarm % dt * dsc_o
          Part % y_o = yi - Part % v * Swarm % dt * dsc_o
          Part % z_o = zi - Part % w * Swarm % dt * dsc_o

          ! Increasing the number of particle reflections
          Swarm % n_reflected(c2) = Swarm % n_reflected(c2) + 1

          ! Debug: write(This_Proc()*1000+k, '(a,2i7, 5(e12.4))')       &
          ! Debug:                         'time, part bounced at ',  &
          ! Debug:                         n, k,                      &
          ! Debug:                         xi, yi, zi, dsc_n, dsc_o

        end if  ! reflected velocity

      end if  ! boundary cell is a wall cell

      !------------------------------------!
      !   The boundary cell is an outlet   !
      !------------------------------------!
      if(Grid % Bnd_Cond_Type(c2) == OUTFLOW  .or.  &
         Grid % Bnd_Cond_Type(c2) == PRESSURE .or.  &
         Grid % Bnd_Cond_Type(c2) == CONVECT) then
        escaped = .true.
        Swarm % n_escaped(c2) = Swarm % n_escaped(c2) + 1
      end if  ! it is an outflow

    end if  ! crossed a boundary cell

  end if  ! s .ne. 0

  end subroutine
