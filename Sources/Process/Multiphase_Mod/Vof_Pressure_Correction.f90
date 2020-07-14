!==============================================================================!
  subroutine Multiphase_Mod_Vof_Pressure_Correction(mult, sol, ini, mass_err)
!------------------------------------------------------------------------------!
!   Correct fluxes on pressure equation due to surface tension and gravity     !
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod,       only: curr_colour => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer                       :: ini
  real                          :: mass_err
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Face_Type),  pointer :: m_flux
  type(Var_Type),   pointer :: vof
  type(Var_Type),   pointer :: u, v, w
  type(Matrix_Type),pointer :: a
  real, contiguous, pointer :: b(:)
  real,             pointer :: u_relax, dt_corr
  integer                   :: c, c1, c2, s
  real                      :: a12, a1_in, a2_in
  real                      :: u_fo, v_fo, w_fo
  real                      :: stens_source, gravity_source, dotprod
  real                      :: dens_f, dens_weight1, dens_weight2, fs
  real                      :: factor1, factor2, correction
!==============================================================================!

  ! Take aliases
  grid      => mult % pnt_grid
  flow      => mult % pnt_flow
  vof       => mult % vof
  u_relax   => flow % u_rel_corr
  dt_corr   => flow % dt_corr
  m_flux    => flow % m_flux
  a         => sol % a
  b         => sol % b % val

  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! Correct for Surface tension
  if(mult % surface_tension > TINY) then

    if (mult % d_func) then
      curr_colour = vof % n
      call Field_Mod_Grad_Variable(flow, vof)
    else
      curr_colour = vof % n
      call Field_Mod_Grad_Variable(flow, vof)
    end if

    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      ! Interpolate VOF gradients
      dotprod =  0.5 * (vof % x(c1) + vof % x(c2)) * grid % dx(s)  &
               + 0.5 * (vof % y(c1) + vof % y(c2)) * grid % dy(s)  &
               + 0.5 * (vof % z(c1) + vof % z(c2)) * grid % dz(s)

      factor1 = u_relax * 0.5 * ( grid % vol(c1) / a % sav(c1)     &
                                + grid % vol(c2) / a % sav(c2) )

      a12 = 0.5 * (mult % curv(c1) + mult % curv(c2)) * mult % surface_tension

      a12 = a12 * factor1 * a % fc(s)

      stens_source = a12 * ( curr_colour(c2) -  curr_colour(c1) - dotprod )

      m_flux % n(s) = m_flux % n(s) + stens_source

      b(c1) = b(c1) - stens_source
      b(c2) = b(c2) + stens_source

    end do

  end if

  ! Correct for Gravity
  if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) >= TINY) then

    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      ! Interpolate gradients
      dotprod =  0.5 * (mult % body_fx(c1) / grid % vol(c1)       &
                      + mult % body_fx(c2) / grid % vol(c2))      &
                      * grid % dx(s)                              &
               + 0.5 * (mult % body_fy(c1) / grid % vol(c1)       &
                      + mult % body_fy(c2) / grid % vol(c2))      &
                      * grid % dy(s)                              &
               + 0.5 * (mult % body_fz(c1) / grid % vol(c1)       &
                      + mult % body_fz(c2) / grid % vol(c2))      &
                      * grid % dz(s)

      gravity_source = ((grid % xf(s) - grid % xc(c1)) * grav_x   &
                      + (grid % yf(s) - grid % yc(c1)) * grav_y   &
                      + (grid % zf(s) - grid % zc(c1)) * grav_z)  &
                      * flow % density(c1)                        &
                      -((grid % xf(s) - grid % xc(c2)) * grav_x   &
                      + (grid % yf(s) - grid % yc(c2)) * grav_y   &
                      + (grid % zf(s) - grid % zc(c2)) * grav_z)  &
                      * flow % density(c2)

      factor1 = u_relax * 0.5 * ( grid % vol(c1) / a % sav(c1)     &
                                + grid % vol(c2) / a % sav(c2) )

      gravity_source =  factor1 * a % fc(s) * (gravity_source - dotprod)

      m_flux % n(s) = m_flux % n(s) + gravity_source

      b(c1) = b(c1) - gravity_source
      b(c2) = b(c2) + gravity_source

    end do

  end if

  ! Introduce temporal correction and subrelaxation
  if (flow % temp_corr) then
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)
      u_fo = fs * u % o(c1) + (1.0 - fs) * u % o(c2)
      v_fo = fs * v % o(c1) + (1.0 - fs) * v % o(c2)
      w_fo = fs * w % o(c1) + (1.0 - fs) * w % o(c2)

      factor2 = u_relax * 0.5 * ( grid % vol(c1) * flow % density(c1)    &
                                  / (a % sav(c1) * dt_corr)              &
                                + grid % vol(c2) * flow % density(c2)    &
                                  / (a % sav(c2) * dt_corr) )

      correction = (1.0 - u_relax) * ( m_flux % star(s) - m_flux % avg(s) )   &
                 + factor2 * ( m_flux % o(s) - ( u_fo * grid % sx(s)          &
                                               + v_fo * grid % sy(s)          &
                                               + w_fo * grid % sz(s) ) )

      m_flux % n(s) = m_flux % n(s) + correction

      b(c1) = b(c1) - correction
      b(c2) = b(c2) + correction

    end do
  end if

  if (mult % phase_change) then
    do c = 1, grid % n_cells
      b(c) = b(c) + mult % flux_rate(c) * grid % vol(c)                    &
                                        * ( 1.0 / mult % phase_dens(1)     &
                                          - 1.0 / mult % phase_dens(2) )
    end do
  end if

  end subroutine
