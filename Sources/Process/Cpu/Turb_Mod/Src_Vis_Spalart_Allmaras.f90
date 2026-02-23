!==============================================================================!
  subroutine Src_Vis_Spalart_Allmaras(Turb, Sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in vis transport equation.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: vis
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: c
  real                       :: x_rat, f_v1, f_v2, f_w, ss, f_t2
  real                       :: dist_v, prod_v, r, gg, dif, dist, chi
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  vis  => Turb % vis
  call Sol % Alias_Native(A, b)

  do c = Cells_In_Domain_And_Buffers()

    if(Turb % model .eq. SPALART_ALLMARAS) then
      dist = Grid % wall_dist(c)
    else if(Turb % model .eq. DES_SPALART) then
      dist = min(Grid % wall_dist(c), Turb % C_des * Turb % h_max(c))
    end if

    !---------------------------------!
    !   Compute the production term   !
    !---------------------------------!
    x_rat = vis % n(c) / (Flow % viscosity(c) / Flow % density(c))
    f_v1  = x_rat**3 / (x_rat**3 + Turb % c_v1**3)
    f_v2  = 1.0 - x_rat/(1.0 + x_rat*f_v1)

    chi   = x_rat

    f_t2  = Turb % c_t3 * exp( -Turb % c_t4 * chi*chi )

    ss    = Flow % vort(c)   &
          + vis % n(c) * f_v2 / (Turb % kappa**2 * dist**2)

    prod_v = Turb % c_b1 * (1.0 - f_t2) * Flow % density(c) * ss * vis % n(c)

    b(c)   = b(c) + prod_v * Grid % vol(c)

    !----------------------------------!
    !   Compute the destruction term   !
    !----------------------------------!
    r  = vis % n(c) / (ss * Turb % kappa**2 * dist**2)
    r  = min(r, 10.0)

    gg = r + Turb % c_w2*(r**6 - r)

    f_w = gg * ((1.0 + Turb % c_w3**6)  &
          / (gg**6 + Turb % c_w3**6))**ONE_SIXTH

    dist_v = ( Turb % c_w1 * f_w                     &
           - Turb % c_b1 / Turb % kappa**2 * f_t2 )  &
           * Flow % density(c) * (vis % n(c) / dist**2)

    A % val(A % dia(c)) = A % val(A % dia(c)) + dist_v * Grid % vol(c)

    !--------------------------------------------!
    !   Compute the first-order diffusion term   !
    !--------------------------------------------!
    dif   = Turb % c_b2                                      &
          * Flow % density(c)                                &
          * (vis % x(c)**2 + vis % y(c)**2 + vis % z(c)**2)  &
          / vis % sigma
    b(c)  = b(c) + dif * Grid % vol(c)

  end do

  end subroutine
