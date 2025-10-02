!==============================================================================!
  subroutine Move_Trapped(Swarm, k)
!------------------------------------------------------------------------------!
!>  Manages the movement of trapped (at an interface) particles within the
!>  swarm. This subroutine is specifically designed for particles that have
!>  interacted with a surface within the computational domain, particularly
!>  for particles adhering to surfaces or those at the interface of different
!>  fluid phases.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Particle-surface interaction: Adjusts the position of trapped particles  !
!     to align with surface norms, ensuring realistic behavior at interfaces.  !
!   * Normal direction movement: Calculates the normal direction to the        !
!     surface at the particle's position and moves the particle accordingly.   !
!   * Position update: Updates the position of trapped particles based on      !
!     the computed movement direction to keep them on the surface.             !
!   * Surface normals: Utilizes the smoothed volume of fluid (VOF) function to !
!     compute surface normals, contributing to the accuracy of the position    !
!     adjustment.                                                              !
!                                                                              !
!   Note                                                                       !
!                                                                              !
!     * This was coppied from Surf % Advance and Surf % Smooth, which may      !
!       sound bad (code duplication), but I had a very strong feeling that     !
!       inheriting the procedure would make things more cumbersome, more       !
!       complex and harder to follow.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm of particles
  integer, intent(in)       :: k      !! particle number (rank)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Field_Type),    pointer :: Flow
  type(Particle_Type), pointer :: Part
  type(Var_Type),      pointer :: u, v, w, smooth
  integer                      :: c
  real                         :: max_dis, rx, ry, rz
  real                         :: nx, ny, nz, dm, smooth_m, smooth_p
! real                         :: r, u_p, v_p, w_p
!==============================================================================!

  ! Take aliases
  Grid   => Swarm % pnt_grid
  Flow   => Swarm % pnt_flow
  Part   => Swarm % Particle(k)
  u      => Flow % u
  v      => Flow % v
  w      => Flow % w
  smooth => Swarm % pnt_vof % smooth
  c      =  Part % cell

  max_dis = 0  ! was used for checking

  !-------------------------------------------------------------------!
  !   Move particles with the velocity of the front.  This strategy   !
  !   outlined here clusters the particles at the trailing end of     !
  !   a bubble.  To make matters worse, it also sheds the particles   !
  !   at sharp(ish) borders of a bubble surface (imagine its skirt)   !
  !-------------------------------------------------------------------!

  ! ! Vector from cell center to old position of the particle
  ! ! (New position was updadet within sub time stepping)
  ! rx = Part % x_o - Grid % xc(c)
  ! ry = Part % y_o - Grid % yc(c)
  ! rz = Part % z_o - Grid % zc(c)
  !
  ! r = sqrt(rx**2 + ry**2 + rz**2)
  ! max_dis = max(max_dis, r)
  !
  ! ! Velocity at the particle
  ! u_p = u % n(c) + u % x(c) * rx  &
  !                + u % y(c) * ry  &
  !                + u % z(c) * rz
  !
  ! v_p = v % n(c) + v % x(c) * rx  &
  !                + v % y(c) * ry  &
  !                + v % z(c) * rz
  !
  ! w_p = w % n(c) + w % x(c) * rx  &
  !                + w % y(c) * ry  &
  !                + w % z(c) * rz
  !
  ! ! New vertex position
  ! Part % x_n = Part % x_o + u_p * Flow % dt
  ! Part % y_n = Part % y_o + v_p * Flow % dt
  ! Part % z_n = Part % z_o + w_p * Flow % dt

  !----------------------------------------------------------------------!
  !   Place them exactly on the surface.  This procedure too, like the   !
  !     one above, clusters particles at a trailing edge of a bubble     !
  !----------------------------------------------------------------------!

  ! Surface normal computed from smoothed vof function
  smooth_m = sqrt(smooth % x(c)**2 + smooth % y(c)**2 + smooth % z(c)**2)
  nx = smooth % x(c) / smooth_m
  ny = smooth % y(c) / smooth_m
  nz = smooth % z(c) / smooth_m

  ! Value at current particle position, updated above
  rx = Part % x_n - Grid % xc(c)
  ry = Part % y_n - Grid % yc(c)
  rz = Part % z_n - Grid % zc(c)
  smooth_p = smooth % n(c) + rx * smooth % x(c)  &
                           + ry * smooth % y(c)  &
                           + rz * smooth % z(c)

  dm = (0.5 - smooth_p)       &
     / (  smooth % x(c)*nx    &
        + smooth % y(c)*ny    &
        + smooth % z(c)*nz)

  rx = dm * nx
  ry = dm * ny
  rz = dm * nz

  ! Move vertex in the surface normal direction
  Part % x_n = Part % x_n + rx
  Part % y_n = Part % y_n + ry
  Part % z_n = Part % z_n + rz

  ! Checked:  rx = Part % x_n - Grid % xc(c)
  ! Checked:  ry = Part % y_n - Grid % yc(c)
  ! Checked:  rz = Part % z_n - Grid % zc(c)
  ! Checked:  smooth_p = smooth % n(c) + rx * smooth % x(c)  &
  ! Checked:                           + ry * smooth % y(c)  &
  ! Checked:                           + rz * smooth % z(c)

  end subroutine
