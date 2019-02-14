!==============================================================================!
  subroutine Numerics_Mod_Advection_Scheme(phi_f,                &
                                           s,                    &
                                           phi,                  &
                                           phi_i, phi_j, phi_k,  &
                                           di, dj, dk,           &
                                           flux)
!------------------------------------------------------------------------------!
!   Computes the value at the cell face using different convective  schemes.   !
!   In this subroutine I try to follow the nomenclature from Basara's and      !
!   Przulj's AIAA paper.                                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real           :: phi_f, phi_f_c, phi_f_u
  integer        :: s
  type(Var_Type) :: phi
  real           :: phi_i(-phi % pnt_grid % n_bnd_cells:  &
                           phi % pnt_grid % n_cells),     &
                    phi_j(-phi % pnt_grid % n_bnd_cells:  &
                           phi % pnt_grid % n_cells),     &
                    phi_k(-phi % pnt_grid % n_bnd_cells:  &
                           phi % pnt_grid % n_cells)
  real           :: di(phi % pnt_grid % n_faces),  &
                    dj(phi % pnt_grid % n_faces),  &
                    dk(phi % pnt_grid % n_faces)
  real           :: flux(phi % pnt_grid % n_faces)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
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
  grid => phi % pnt_grid

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
    phi_star = phi % n(d) - 2.0 * (  phi_i(c)*di(s)   &
                                   + phi_j(c)*dj(s)   &
                                   + phi_k(c)*dk(s) )
  else
    phi_star = phi % n(d) + 2.0 * (  phi_i(c)*di(s)   &
                                   + phi_j(c)*dj(s)   &
                                   + phi_k(c)*dk(s) )
  end if

  phi_u = max( phi % min(c), min(phi_star, phi % max(c)) )

  rj = ( phi % n(c) - phi_u ) / ( phi % n(d)-phi % n(c) + 1.0e-16 )

  g_d = 0.5 * fj * (1.0+fj)
  g_u = 0.5 * fj * (1.0-fj)

  if(phi % adv_scheme .eq. CENTRAL) then
    phij = fj

  else if(phi % adv_scheme .eq. QUICK) then
    rj = ( phi % n(c) - phi_u ) / ( phi % n(d)-phi % n(c) + 1.0e-12 )
    alfa = 0.0
    phij = (g_d - alfa) + (g_u + alfa) * rj

  else if(phi % adv_scheme .eq. LUDS) then
    alfa = 0.5 * fj * (1+fj)
    phij = (g_d - alfa) + (g_u + alfa) * rj

  else if(phi % adv_scheme .eq. MINMOD) then
    phij = fj * max(0.0, min(rj,1.0))

  else if(phi % adv_scheme .eq. SMART) then
    beta1 = 3.0
    beta2 = 1.0
    phij = max( 0.0, min( (beta1-1.0)*rj, g_d+g_u*rj, beta2 ) )

  else if(phi % adv_scheme .eq. AVL_SMART) then
    beta1 = 1.0 + fj*(2.0+fj) 
    beta2 = fj*(2.0-fj) 
    phij = max( 0.0, min( (beta1-1.0)*rj, g_d+g_u*rj, beta2 ) )

  else if(phi % adv_scheme .eq. SUPERBEE) then
    phij = 0.5 * max( 0.0, min( 2.0*rj,1.0 ), min( rj,2.0 ) )

  else if(phi % adv_scheme .eq. UPWIND) then
    phi_f = phi % n(c)  ! upwind value
    return

  else if(phi % adv_scheme .eq. BLENDED) then
    phi_f_c = phi % n(c) + fj * sign * (phi % n(c2)-phi % n(c1))  ! central part
    phi_f_u = phi % n(c)                                          ! upwind part
    ! Blended value
    phi_f   =        phi % blend  * phi_f_c   &
            + (1.0 - phi % blend) * phi_f_u
    return
  end if

  phi_f = phi % n(c) + phij * sign * (phi % n(c2)-phi % n(c1))

  if(phi % adv_scheme .eq. GAMMA) then
    beta = 0.1

    if(flux(s) > 0.0) then
      phi_star = 1.0 - (phi % n(d) - phi % n(c))/(2.0 * (  phi_i(c)*di(s)   &
                                                         + phi_j(c)*dj(s)   &
                                                         + phi_k(c)*dk(s)))
    else
      phi_star = 1.0 + (phi % n(d) - phi % n(c))/(2.0 * (  phi_i(c)*di(s)   &
                                                         + phi_j(c)*dj(s)   &
                                                         + phi_k(c)*dk(s)))
    end if

    gamma_c = phi_star / beta

    if(phi_star < beta.and.phi_star > 0.0) then
      phi_f = (1.0 - gamma_c*(1.0 - grid % f(s))) * phi % n(c)   &
                   + gamma_c*(1.0 - grid % f(s))  * phi % n(d)
    else if(phi_star < 1.0.and.phi_star >= beta) then
       phi_f =        grid % f(s)  * phi % n(c)  &
             + (1.0 - grid % f(s)) * phi % n(d)
    else
      phi_f = phi % n(c)
    end if
  end if

  end subroutine
