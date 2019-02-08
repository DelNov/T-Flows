!==============================================================================!
  subroutine Numerics_Mod_Advection_Scheme(flow,                             &
                                           phi_f, s,                         &
                                           phi, phi_min, phi_max,            &
                                           phi_i, phi_j, phi_k, di, dj, dk,  &
                                           scheme, blend)
!------------------------------------------------------------------------------!
!   Computes the value at the cell face using different convective  schemes.   !
!   In this subroutine I try to follow the nomenclature from Basara's and      !
!   Przulj's AIAA paper.                                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod,  only: Field_Type
  use Grid_Mod,   only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  real                     :: phi_f, phi_f_c, phi_f_u
  integer                  :: s
  real                     :: phi    (-flow % pnt_grid % n_bnd_cells:  &
                                       flow % pnt_grid % n_cells),     &
                              phi_min(-flow % pnt_grid % n_bnd_cells:  &
                                       flow % pnt_grid % n_cells),     &
                              phi_max(-flow % pnt_grid % n_bnd_cells:  &
                                       flow % pnt_grid % n_cells)
  real                     :: phi_i(-flow % pnt_grid % n_bnd_cells:  &
                                     flow % pnt_grid % n_cells),     &
                              phi_j(-flow % pnt_grid % n_bnd_cells:  &
                                     flow % pnt_grid % n_cells),     &
                              phi_k(-flow % pnt_grid % n_bnd_cells:  &
                                     flow % pnt_grid % n_cells)
  real                     :: di(flow % pnt_grid % n_faces),  &
                              dj(flow % pnt_grid % n_faces),  &
                              dk(flow % pnt_grid % n_faces)
  integer                  :: scheme  ! advection scheme
  real                     :: blend
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  real,            pointer :: flux(:)
  integer                  :: c1, c2, c, d
  real                     :: fj ! flow oriented interpolation factor
  real                     :: g_d, g_u, alfa, beta1, beta2 
  real                     :: phij, phi_u, phi_star, rj, sign, gamma_c, beta
!==============================================================================!
!
!               Flux > 0
!   +---------+---------+---------+---------+
!   |         |         |         |         |
!   |         |   c1    |   c2    |         |
!   |    o    |    o  ==s=>  o    |    o    |----> xi
!   |    U    |    c    |    d    |         |
!   |         |         |         |         |
!   +---------+---------+---------+---------+   
!
!
!               Flux < 0
!   +---------+---------+---------+---------+
!   |         |         |         |         |
!   |         |   c1    |   c2    |         |
!   |    o    |    o  <=s==  o    |    o    |----> xi
!   |         |    d    |    c    |    U    |
!   |         |         |         |         |
!   +---------+---------+---------+---------+   
!
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  flux => flow % flux

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)

  if(flux(s) > 0.0) then ! goes from c1 to c2
    fj   = 1.0 - grid % f(s)
    c    = c1
    d    = c2
    sign = +1.0
  else ! flux(s) < 0.0   ! goes from c2 to c1
    fj = grid % f(s)
    c    = c2
    d    = c1
    sign = -1.0
  end if

  if(flux(s) > 0.0) then
    phi_star = phi(d) - 2.0 * ( phi_i(c)*di(s) &
                               +phi_j(c)*dj(s) &
                               +phi_k(c)*dk(s) )
  else
    phi_star = phi(d) + 2.0 * ( phi_i(c)*di(s) &
                               +phi_j(c)*dj(s) &
                               +phi_k(c)*dk(s) )
  end if

  phi_u = max( phi_min(c), min(phi_star, phi_max(c)) )

  rj = ( phi(c) - phi_u ) / ( phi(d)-phi(c) + 1.0e-16 )

  g_d = 0.5 * fj * (1.0+fj)
  g_u = 0.5 * fj * (1.0-fj)

  if(scheme .eq. CENTRAL) then
    phij = fj

  else if(scheme .eq. QUICK) then
    rj = ( phi(c) - phi_u ) / ( phi(d)-phi(c) + 1.0e-12 )
    alfa = 0.0
    phij = (g_d - alfa) + (g_u + alfa) * rj

  else if(scheme .eq. LUDS) then
    alfa = 0.5 * fj * (1+fj)
    phij = (g_d - alfa) + (g_u + alfa) * rj

  else if(scheme .eq. MINMOD) then
    phij = fj * max(0.0, min(rj,1.0))

  else if(scheme .eq. SMART) then
    beta1 = 3.0
    beta2 = 1.0
    phij = max( 0.0, min( (beta1-1.0)*rj, g_d+g_u*rj, beta2 ) )

  else if(scheme .eq. AVL_SMART) then
    beta1 = 1.0 + fj*(2.0+fj) 
    beta2 = fj*(2.0-fj) 
    phij = max( 0.0, min( (beta1-1.0)*rj, g_d+g_u*rj, beta2 ) )

  else if(scheme .eq. SUPERBEE) then
    phij = 0.5 * max( 0.0, min( 2.0*rj,1.0 ), min( rj,2.0 ) )

  else if(scheme .eq. UPWIND) then
    phi_f = phi(c)  ! upwind value
    return

  else if(scheme .eq. BLENDED) then
    phi_f_c = phi(c) + fj * sign * (phi(c2)-phi(c1))  ! central part
    phi_f_u = phi(c)                                  ! upwind part
    phi_f = blend * phi_f_c + (1.0-blend) * phi_f_u   ! blended value
    return
  end if

  phi_f = phi(c) + phij * sign * (phi(c2)-phi(c1))

  if(scheme .eq. GAMMA) then
    beta = 0.1

    if(flux(s) > 0.0) then
      phi_star = 1.0 - (phi(d) - phi(c))/(2.0 * ( phi_i(c)*di(s) &
                                                + phi_j(c)*dj(s) &
                                                + phi_k(c)*dk(s)))
    else
      phi_star = 1.0 + (phi(d) - phi(c))/(2.0 * ( phi_i(c)*di(s) &
                                                + phi_j(c)*dj(s) &
                                                + phi_k(c)*dk(s)))
    end if

    gamma_c = phi_star / beta

    if(phi_star < beta.and.phi_star > 0.0) then
      phi_f = (1.0 - gamma_c*(1.0 - grid % f(s))) * phi(c)   &
                   + gamma_c*(1.0 - grid % f(s))  * phi(d)
    else if(phi_star < 1.0.and.phi_star >= beta) then
       phi_f =        grid % f(s)  * phi(c)  &
             + (1.0 - grid % f(s)) * phi(d)
    else
      phi_f = phi(c)
    end if
  end if

  end subroutine
