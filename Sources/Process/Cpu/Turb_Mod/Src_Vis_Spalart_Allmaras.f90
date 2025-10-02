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
  real                       :: x_rat, f_v1, f_v2, f_w, ss
  real                       :: dist_v, prod_v, r, gg, dif, dist
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  vis  => Turb % vis
  call Sol % Alias_Native(A, b)

  if(Turb % model .eq. SPALART_ALLMARAS) then

    do c = Cells_In_Domain_And_Buffers()

      dist = Grid % wall_dist(c)

      !---------------------------------!
      !   Compute the production term   !
      !---------------------------------!
      x_rat  = vis % n(c) / (Flow % viscosity(c) / Flow % density(c))
      f_v1   = x_rat**3 / (x_rat**3 + Turb % c_v1**3)
      f_v2   = 1.0 - x_rat/(1.0 + x_rat*f_v1)
      ss     = Flow % vort(c)   &
             + vis % n(c) * f_v2 / (Turb % kappa**2 * dist**2)
      prod_v = Turb % c_b1 * Flow % density(c) * ss * vis % n(c)
      b(c)   = b(c) + prod_v * Grid % vol(c)

      !----------------------------------!
      !   Compute the destruction term   !
      !----------------------------------!
      r      = vis % n(c) / (ss * Turb % kappa**2 * dist**2)
      gg     = r + Turb % c_w2*(r**6 - r)
      f_w    = gg*((1.0 + Turb % c_w3**6)  &
             / (gg**6 + Turb % c_w3**6))**ONE_SIXTH
      dist_v = Turb % c_w1 * Flow % density(c) * f_w  &
             * (vis % n(c) / dist**2)
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

  else if(Turb % model .eq. DES_SPALART) then
    do c = Cells_In_Domain_And_Buffers()

      ! What is 0.65 here?  A ghost number
      dist = min(Grid % wall_dist(c), 0.65 * Turb % h_max(c))

      !---------------------------------!
      !   Compute the production term   !
      !---------------------------------!
      x_rat  = vis % n(c) / (Flow % viscosity(c) / Flow % density(c))
      f_v1   = x_rat**3 / (x_rat**3 + Turb % c_v1**3)
      f_v2   = 1.0 - x_rat/(1.0 + x_rat*f_v1)
      ss     = Flow % vort(c) + vis % n(c) * f_v2 / (Turb % kappa**2 * dist**2)
      prod_v = Turb % c_b1 * Flow % density(c) * ss * vis % n(c)
      b(c)   = b(c) + prod_v * Grid % vol(c)

      !-----------------------------------!
      !   Compute the destruction  term   !
      !-----------------------------------!
      r      = vis % n(c) / (ss * Turb % kappa**2 * dist**2)
      gg     = r + Turb % c_w2 * (r**6 - r)
      f_w    = gg*((1.0 + Turb % c_w3**6)  &
             / (gg**6 + Turb % c_w3**6))**ONE_SIXTH
      dist_v = Turb % c_w1 * Flow % density(c) * f_w * (vis % n(c) / dist**2)
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
  end if

  end subroutine
