!==============================================================================!
  subroutine Multiphase_Mod_Vof_Pressure_Correction(mult, sol)
!------------------------------------------------------------------------------!
!   Correct fluxes on pressure equation due to surface tension and gravity     !
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod,       only: dens_factor => r_face_01,   &
                            curr_colour => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Face_Type),  pointer :: m_flux
  type(Var_Type),   pointer :: vof
  type(Var_Type),   pointer :: dist_func
  type(Matrix_Type),pointer :: a
  real,             pointer :: b(:)
  integer                   :: c1, c2, s
  real                      :: a12, a1_in, a2_in, df, cf, dt
  real                      :: stens_source, gravity_source, dotprod
!==============================================================================!

  ! Take aliases
  grid      => mult % pnt_grid
  flow      => mult % pnt_flow
  m_flux    => flow % m_flux
  vof       => mult % vof
  dist_func => mult % dist_func
  a         => sol % a
  b         => sol % b % val

  dt = flow % dt

  if(mult % surface_tension > TINY) then

    if (mult % d_func) then
      curr_colour = dist_func % oo
    else
      curr_colour = vof % n
    end if

    if (mult % n_conv_norm > 0) then
      call Field_Mod_Grad_Variable(flow, vof)
    end if

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Face is inside the domain
      if(c2 > 0) then

        a1_in = flow % density(c1) * grid % vol(c1) / dt
        a2_in = flow % density(c2) * grid % vol(c2) / dt

        df = 0.5 * ( grid % vol(c1) / (a % sav(c1) - a1_in)    &
                   + grid % vol(c2) / (a % sav(c2) - a2_in) )

        cf = flow % density_f(s) / dt

        ! Interpolate VOF gradients
        dotprod =  0.5 * (vof % x(c1) + vof % x(c2)) * grid % dx(s)  &
                 + 0.5 * (vof % y(c1) + vof % y(c2)) * grid % dy(s)  &
                 + 0.5 * (vof % z(c1) + vof % z(c2)) * grid % dz(s)

!        a12 = 0.5 * a % fc(s)                    &
!             * ( grid % vol(c1) / a % sav(c1)    &
!               + grid % vol(c2) / a % sav(c2) )
!
!        a12 = a12 * 0.5 * ( mult % curv(c1) + mult % curv(c2) )  &
!                          * mult % surface_tension
        a12 = 0.5 * (mult % curv(c1) + mult % curv(c2)) * mult % surface_tension

        a12 = a12 * df * a % fc(s) / (1.0 + cf * df)

        stens_source = a12 * ( curr_colour(c2) -  curr_colour(c1) - dotprod )

        m_flux % n(s) = m_flux % n(s) + stens_source

        b(c1) = b(c1) - stens_source
        b(c2) = b(c2) + stens_source

      end if

    end do

  end if

  if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) >= TINY) then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Face is inside the domain
      if(c2 > 0) then

        ! Interpolate density 
        a1_in = flow % density(c1) *  grid % vol(c1) / dt
        a2_in = flow % density(c2) *  grid % vol(c2) / dt

        df = 0.5 * ( grid % vol(c1) / (a % sav(c1) - a1_in)    &
                   + grid % vol(c2) / (a % sav(c2) - a2_in) )

        cf = flow % density_f(s) / dt

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

!        a12 = 0.5 * a % fc(s)                    &
!             * ( grid % vol(c1) / a % sav(c1)    &
!               + grid % vol(c2) / a % sav(c2) )

!        gravity_source =  a12 * (gravity_source - dotprod)
        gravity_source =  df * a % fc(s) / (1.0 + cf * df)  &
                        * (gravity_source - dotprod)

        m_flux % n(s) = m_flux % n(s) + gravity_source

        b(c1) = b(c1) - gravity_source
        b(c2) = b(c2) + gravity_source

      end if

    end do

  end if

  end subroutine
