!==============================================================================!
  subroutine Particles_Mod_Interpolate_Velocity(flow, part)
!------------------------------------------------------------------------------!
!   Interpolates velocity at the particle.                                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod, only: Field_Type
  use Grid_Mod,  only: Grid_Type
  use Var_Mod,   only: Var_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: flow
  type(Particle_Type)         :: part
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: s, c1, c2, c
  real                     :: rx, ry, rz  ! vector connecting paticle and cell
  real                     :: xc, yc, zc  ! cell center coordinates
  real                     :: up, vp, wp  ! velocity at particle position
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w

  !-------------------------------------------------!
  !   Browse through all cells as another example   !
  !-------------------------------------------------!
  do c = 1, grid % n_cells

    ! Cell cenre coorindates
    xc = grid % xc(c)
    yc = grid % yc(c)
    zc = grid % zc(c)

    ! Vector which connects particle position and cell centre
    rx = part % x - xc
    ry = part % y - yc
    rz = part % z - zc

    ! Compute velocities at the particle from velocity gradients
    up = u % n(c)       &  ! u velocity at the new time step (% n)
       + u % x(c) * rx  &  ! u % x is gradient du/dx
       + u % y(c) * ry  &  ! u % y is gradient du/dy
       + u % z(c) * rz     ! u % x is gradient du/dz

    vp = v % n(c)       &  ! v velocity at the new time step (% n)
       + v % x(c) * rx  &  ! v % x is gradient dv/dx
       + v % y(c) * ry  &  ! v % y is gradient dv/dy
       + v % z(c) * rz     ! v % x is gradient dv/dz

    wp = w % n(c)       &  ! w velocity at the new time step (% n)
       + w % x(c) * rx  &  ! w % x is gradient dw/dx
       + w % y(c) * ry  &  ! w % y is gradient dw/dy
       + w % z(c) * rz     ! w % x is gradient dw/dz

  end do

  end subroutine
