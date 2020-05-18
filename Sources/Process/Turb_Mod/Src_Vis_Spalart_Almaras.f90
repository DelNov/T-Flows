!==============================================================================!
  subroutine Turb_Mod_Src_Vis_Spalart_Almaras(turb, sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in vis transport equation.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: vis
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c
  real                       :: x_rat, f_v1, f_v2, f_w, ss
  real                       :: dist_v, prod_v, r, gg, dif, dist
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  vis  => turb % vis
  call Solver_Mod_Alias_System(sol, a, b)

  if(turb % model .eq. SPALART_ALLMARAS) then

    do c = 1, grid % n_cells

      !---------------------------------!
      !   Compute the production term   !
      !---------------------------------!
      x_rat  = vis % n(c) / flow % viscosity(c)
      f_v1   = x_rat**3/(x_rat**3 + c_v1**3)
      f_v2   = 1.0 - x_rat/(1.0 + x_rat*f_v1)
      ss     = flow % vort(c)   &
             + vis % n(c) * f_v2 / (kappa**2 * grid % wall_dist(c)**2)
      prod_v = c_b1 * flow % density(c) * ss * vis % n(c)
      b(c)   = b(c) + prod_v * grid % vol(c)

      !----------------------------------!
      !   Compute the destruction term   !
      !----------------------------------!
      r      = vis % n(c)/(ss * kappa**2 * grid % wall_dist(c)**2)
      gg     = r + c_w2*(r**6 - r)
      f_w    = gg*((1.0 + c_w3**6)/(gg**6 + c_w3**6))**(1.0/6.0)
      dist_v = c_w1 * flow % density(c) * f_w  &
             * (vis % n(c) / grid % wall_dist(c)**2)
      a % val(a % dia(c)) = a % val(a % dia(c)) + dist_v * grid % vol(c)

      !--------------------------------------------!
      !   Compute the first-order diffusion term   !
      !--------------------------------------------!
      dif   = c_b2                                       &
            * flow % density(c)                          &
            * (vis % x(c) + vis % y(c) + vis % z(c))**2  &
            / vis % sigma
      b(c)  = b(c) + dif * grid % vol(c)
    end do

  else if(turb % model .eq. DES_SPALART) then
    do c = 1, grid % n_cells

      ! What is 0.65 here?  A ghost number
      dist = min(grid % wall_dist(c), 0.65 * turb % h_max(c))

      !---------------------------------!
      !   Compute the production term   !
      !---------------------------------!
      x_rat  = vis % n(c) / flow % viscosity(c)
      f_v1   = x_rat**3 / (x_rat**3 + c_v1**3)
      f_v2   = 1.0 - x_rat/(1.0 + x_rat*f_v1)
      ss     = flow % vort(c) + vis % n(c) * f_v2 / (kappa**2 * dist**2)
      prod_v = c_b1 * flow % density(c) * ss * vis % n(c)
      b(c)   = b(c) + prod_v * grid % vol(c)

      !-----------------------------------!
      !   Compute the destruction  term   !
      !-----------------------------------!
      r      = vis % n(c) / (ss * kappa**2 * dist**2)
      gg     = r + c_w2 * (r**6 - r)
      f_w    = gg*((1.0 + c_w3**6) / (gg**6 + c_w3**6))**(1.0/6.0)
      dist_v = c_w1 * flow % density(c) * f_w * (vis % n(c) / dist**2)
      a % val(a % dia(c)) = a % val(a % dia(c)) + dist_v * grid % vol(c)

      !--------------------------------------------!
      !   Compute the first-order diffusion term   !
      !--------------------------------------------!
      dif   = c_b2                                       &
            * flow % density(c)                          &
            * (vis % x(c) + vis % y(c) + vis % z(c))**2  &
            / vis % sigma
      b(c)  = b(c) + dif * grid % vol(c)
    end do
  end if

  end subroutine
