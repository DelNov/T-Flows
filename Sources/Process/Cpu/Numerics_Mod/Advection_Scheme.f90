!==============================================================================!
  subroutine Numerics_Mod_Advection_Scheme(phi_f,                &
                                           s,                    &
                                           phi,                  &
                                           phi_min,              &
                                           phi_max,              &
                                           flux)
!------------------------------------------------------------------------------!
!>  The subroutine Numerics_Mod_Advection_Scheme within the Numerics_Mod module
!>  calculates the value at the cell face using various advection schemes. The
!>  subroutine follows the nomenclature from Basara's and Przulj's AIAA paper,
!>  establishing a systematic approach to implement different advection schemes.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Face value computation: Determines the interpolated value of a           !
!     transported variable at a face, crucial for convective flux calculations.!
!   * Scheme Flexibility: Supports multiple advection schemes like CENTRAL,    !
!     QUICK, LUDS, MINMOD, SMART, AVL_SMART, SUPERBEE, UPWIND, BLENDED, and    !
!     GAMMA, allowing for customization based on flow characteristics and      !
!     accuracy requirements.                                                   !
!   * Flux direction consideration: Accounts for the direction of flux across  !
!     the face, ensuring appropriate upwind or downwind biasing as needed.     !
!   * Interpolation: Employs different interpolation techniques (upwind,       !
!     central, blended, etc.) to compute face values, balancing accuracy       !
!     and stability.                                                           !
!   * Boundedness and limiting: Incorporates limiting functions to maintain    !
!     boundedness of the solution, preventing non-physical oscillations.       !
!   * Gamma scheme support: Specifically implements the Gamma scheme for       !
!     enhanced accuracy in certain flow regimes.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real                       :: phi_f  !! interpolatied value at the face
  integer,        intent(in) :: s      !! face number
  type(Var_Type), intent(in) :: phi    !! transported variable
  real,           intent(in) :: phi_min(-phi % pnt_grid % n_bnd_cells:  &
                                         phi % pnt_grid % n_cells)
                                      !! lower bound array for phi
  real,           intent(in) :: phi_max(-phi % pnt_grid % n_bnd_cells:  &
                                         phi % pnt_grid % n_cells)
                                      !! lower bound array for phi
  real,           intent(in) :: flux(phi % pnt_grid % n_faces)
                                       !! volume flux at faces
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c1, c2, c, d
  real                     :: fj       ! flow oriented interpolation factor
  real                     :: g_d, g_u, alfa, beta1, beta2, denom
  real                     :: phi_f_c, phi_f_u
  real                     :: phij, phi_u, phi_star, rj, sgn, gamma_c, beta
!------------------------------------------------------------------------------!
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
!==============================================================================!

  ! Take aliases
  Grid => phi % pnt_grid

  c1 = Grid % faces_c(1,s)
  c2 = Grid % faces_c(2,s)

  if(flux(s) > 0.0) then ! goes from c1 to c2
    fj  = 1.0 - Grid % f(s)
    c   = c1
    d   = c2
    sgn = +1.0
  else ! flux(s) < 0.0   ! goes from c2 to c1
    fj  = Grid % f(s)
    c   = c2
    d   = c1
    sgn = -1.0
  end if

  if(flux(s) > 0.0) then
    phi_star = phi % n(d) - 2.0 * (  phi % x(c)*Grid % dx(s)   &
                                   + phi % y(c)*Grid % dy(s)   &
                                   + phi % z(c)*Grid % dz(s) )
  else
    phi_star = phi % n(d) + 2.0 * (  phi % x(c)*Grid % dx(s)   &
                                   + phi % y(c)*Grid % dy(s)   &
                                   + phi % z(c)*Grid % dz(s) )
  end if

  phi_u = max( phi_min(c), min(phi_star, phi_max(c)) )

  denom = phi % n(d) - phi % n(c)
  if( abs(denom) < FEMTO ) then
    denom = denom + sign(FEMTO, denom)  ! to avoid changing sign
  endif
  rj = (phi % n(c) - phi_u) / denom

  g_d = 0.5 * fj * (1.0+fj)
  g_u = 0.5 * fj * (1.0-fj)

  if(phi % adv_scheme .eq. CENTRAL) then
    phij = fj

  else if(phi % adv_scheme .eq. QUICK) then
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
    phij = max( 0.0, min( (beta1-1.0)*rj, g_d + g_u*rj, beta2 ) )

  else if(phi % adv_scheme .eq. AVL_SMART) then
    beta1 = 1.0 + fj*(2.0+fj) 
    beta2 = fj*(2.0-fj) 
    phij = max( 0.0, min( (beta1-1.0)*rj, g_d + g_u*rj, beta2 ) )

  else if(phi % adv_scheme .eq. SUPERBEE) then
    phij = 0.5 * max( 0.0, min( 2.0*rj,1.0 ), min( rj,2.0 ) )

  else if(phi % adv_scheme .eq. UPWIND) then
    phi_f = phi % n(c)  ! upwind value
    return

  else if(phi % adv_scheme .eq. BLENDED) then
    phi_f_c = phi % n(c) + fj * sgn * (phi % n(c2)-phi % n(c1))  ! central part
    phi_f_u = phi % n(c)                                         ! upwind part
    ! Blended value
    phi_f   =        phi % blend  * phi_f_c   &
            + (1.0 - phi % blend) * phi_f_u
    return
  end if

  phi_f = phi % n(c) + phij * sgn * (phi % n(c2)-phi % n(c1))

  if(phi % adv_scheme .eq. GAMMA) then
    beta = 0.1

    if(flux(s) > 0.0) then
      phi_star = 1.0 - (phi % n(d) - phi % n(c))              &
                     / (2.0 * (  phi % x(c) * Grid % dx(s)    &
                               + phi % y(c) * Grid % dy(s)    &
                               + phi % z(c) * Grid % dz(s)))
    else
      phi_star = 1.0 + (phi % n(d) - phi % n(c))              &
                     / (2.0 * (  phi % x(c) * Grid % dx(s)    &
                               + phi % y(c) * Grid % dy(s)    &
                               + phi % z(c) * Grid % dz(s)))
    end if

    gamma_c = phi_star / beta

    if(phi_star < beta .and. phi_star > 0.0) then
      phi_f = (1.0 - gamma_c*(1.0 - Grid % f(s))) * phi % n(c)   &
                   + gamma_c*(1.0 - Grid % f(s))  * phi % n(d)
    else if(phi_star < 1.0.and.phi_star >= beta) then
       phi_f =        Grid % f(s)  * phi % n(c)  &
             + (1.0 - Grid % f(s)) * phi % n(d)
    else
      phi_f = phi % n(c)
    end if
  end if

  end subroutine
