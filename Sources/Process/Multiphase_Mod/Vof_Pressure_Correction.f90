!==============================================================================!
  subroutine Multiphase_Mod_Vof_Pressure_Correction(mult, sol)
!------------------------------------------------------------------------------!
!   Correct fluxes on pressure equation due to surface tension and gravity     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Face_Type),   pointer :: v_flux
  type(Var_Type),    pointer :: col
  type(Var_Type),    pointer :: u, v, w
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:)
  real,              pointer :: u_relax, dt_corr
  integer                    :: c, c1, c2, s, nb, nc
  real                       :: a12, fs
  real                       :: u_fo, v_fo, w_fo
  real                       :: stens_source, dotprod
  real                       :: factor2, correction, dens_h, curv_f
!==============================================================================!

  ! Take aliases
  grid    => mult % pnt_grid
  flow    => mult % pnt_flow
  ! col     => mult % smooth
  col     => mult % vof
  u_relax => flow % u_rel_corr
  dt_corr => flow % dt_corr
  v_flux  => flow % v_flux
  a       => sol % a
  b       => sol % b % val

  nb = grid % n_bnd_cells
  nc = grid % n_cells

  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! Correct for Surface tension
  if(mult % surface_tension > TINY) then

    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      ! Interpolate VOF gradients
      dens_h = 2.0 / ( 1.0 / flow % density(c1) + 1.0 / flow % density(c2) )

      ! Unit for dotprod: [1/m]
      dotprod = 0.5 * dens_h                                                 &
                    * ( mult % curv(c1) * col % x(c1) / flow % density(c1)   &
                      + mult % curv(c2) * col % x(c2) / flow % density(c2) ) &
                    * grid % dx(s)                                           &
              + 0.5 * dens_h                                                 &
                    * ( mult % curv(c1) * col % y(c1) / flow % density(c1)   &
                      + mult % curv(c2) * col % y(c2) / flow % density(c2) ) &
                    * grid % dy(s)                                           &
              + 0.5 * dens_h                                                 &
                    * ( mult % curv(c1) * col % z(c1) / flow % density(c1)   &
                      + mult % curv(c2) * col % z(c2) / flow % density(c2) ) &
                    * grid % dz(s)

      ! Unit for a12: [m^4s/kg]
      a12 = u_relax * 0.5 * ( grid % vol(c1) / a % sav(c1)     &
                            + grid % vol(c2) / a % sav(c2) ) * a % fc(s)

      ! Curvature at the face; unit: [1/m]
      curv_f = 0.5 * ( mult % curv(c1) + mult % curv(c2) )

      ! Unit for stens_source: [kg/s^2 * m^4s/kg * 1/m = m^3/s]
      stens_source = mult % surface_tension * a12               &
                   * ( curv_f * (col % n(c2) -  col % n(c1)) - dotprod )

      v_flux % n(s) = v_flux % n(s) + stens_source

      b(c1) = b(c1) - stens_source
      b(c2) = b(c2) + stens_source

    end do

  end if

  ! Introduce temporal correction and subrelaxation
  ! (See equation 3.61 in Denner's thesis)
  if (flow % temp_corr) then
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)
      u_fo = fs * u % o(c1) + (1.0 - fs) * u % o(c2)
      v_fo = fs * v % o(c1) + (1.0 - fs) * v % o(c2)
      w_fo = fs * w % o(c1) + (1.0 - fs) * w % o(c2)

      dens_h = 2.0 / ( 1.0 / flow % density(c1) + 1.0 / flow % density(c2) )

      ! Unit for factor2: [1]
      factor2 = u_relax * 0.5 * dens_h * ( grid % vol(c1)               &
                                         / (a % sav(c1) * dt_corr)      &
                                         + grid % vol(c2)               &
                                         / (a % sav(c2) * dt_corr) )

      ! Unit for correction [m^3/s]
      correction = (1.0 - u_relax) * ( v_flux % star(s) - v_flux % avg(s) )   &
                 + factor2 * ( v_flux % o(s) - ( u_fo * grid % sx(s)          &
                                               + v_fo * grid % sy(s)          &
                                               + w_fo * grid % sz(s) ) )

      v_flux % n(s) = v_flux % n(s) + correction

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
