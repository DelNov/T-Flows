!==============================================================================!
  subroutine Swarm_Mod_Check_Periodicity(swarm, k)
!------------------------------------------------------------------------------!
!   Check if particle left the cell (which means periodic domain as well)      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  real                         :: xn, yn, zn   ! node coordinates
  real                         :: xc, yc, zc   ! cell coordinates
  real                         :: rxc, ryc, rzc, rxn, ryn, rzn
!==============================================================================!

  ! Take aliases
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid
  part => swarm % particle(k)

  if(grid % per_x < TINY .and.  &
     grid % per_y < TINY .and.  &
     grid % per_z < TINY) return

  ! X-direction
  xn  = grid % xn(part % node)
  xc  = grid % xc(part % cell)
  rxn = part % x_n - xn
  rxc = xc - part % x_n
  if(rxc*rxn < 0.0) then
    if(flow % u % n(part % cell) > 0.0) then
      part % x_n = part % x_n - grid % per_x
      part % x_o = part % x_o - grid % per_x
    else
      part % x_n = part % x_n + grid % per_x
      part % x_o = part % x_o + grid % per_x
    end if
  end if

  ! Y-direction
  yn = grid % yn(part % node)
  yc = grid % yc(part % cell)
  ryn = part % y_n - yn
  ryc = yc - part % y_n
  if(ryc*ryn < 0.0) then
    if(flow % v % n(part % cell) > 0.0) then
      part % y_n = part % y_n - grid % per_y
      part % y_o = part % y_o - grid % per_y
    else
      part % y_n = part % y_n + grid % per_y
      part % y_o = part % y_o + grid % per_y
    end if
  end if

  ! Z-direction
  zn = grid % zn(part % node)
  zc = grid % zc(part % cell)
  rzn = part % z_n - zn
  rzc = zc - part % z_n
  if(rzc*rzn < 0.0) then
    if(flow % w % n(part % cell) > 0.0) then
      part % z_n = part % z_n - grid % per_z
      part % z_o = part % z_o - grid % per_z
    else
      part % z_n = part % z_n + grid % per_z
      part % z_o = part % z_o + grid % per_z
    end if
  end if

  if(rxc*rxn < 0.0 .or.  &
     ryc*ryn < 0.0 .or.  &
     rzc*rzn < 0.0) then
    part % node = 0
  end if

  end subroutine
