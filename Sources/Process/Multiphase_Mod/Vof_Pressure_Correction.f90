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
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Face_Type),  pointer :: m_flux
  type(Var_Type),   pointer :: vof, u
  type(Matrix_Type),pointer :: a
  real,             pointer :: b(:)
  integer                   :: c1, c2, s
  real                      :: a12
  real                      :: stens_source, gravity_source, dotprod
!==============================================================================!

  ! Take aliases
  grid   => mult % pnt_grid
  flow   => mult % pnt_flow
  m_flux => flow % m_flux
  u      => flow % u
  vof    => mult % vof
  a      => sol % a
  b      => sol % b % val

  if(abs(surface_tension) > TINY) then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Face is inside the domain
      if(c2 > 0) then

        ! Interpolate VOF gradients
        dotprod =  0.5 * (vof % x(c1) + vof % x(c2)) * grid % dx(s)  &
                 + 0.5 * (vof % y(c1) + vof % y(c2)) * grid % dy(s)  &
                 + 0.5 * (vof % z(c1) + vof % z(c2)) * grid % dz(s)

        a12 = 0.5 * a % fc(s)                    &
             * ( grid % vol(c1) / a % sav(c1)    &
               + grid % vol(c2) / a % sav(c2) )

        a12 = a12 * 0.5 * ( mult % curv(c1) + mult % curv(c2) )  &
                          * surface_tension

        stens_source = a12 * ( vof % n(c2) -  vof % n(c1) - dotprod )

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

        ! Interpolate gradients
        dotprod =  0.5 * (body_fx(c1) / grid % vol(c1)              &
                        + body_fx(c2) / grid % vol(c2))             &
                        * grid % dx(s)                              &
                 + 0.5 * (body_fy(c1) / grid % vol(c1)              &
                        + body_fy(c2) / grid % vol(c2))             &
                        * grid % dy(s)                              &
                 + 0.5 * (body_fz(c1) / grid % vol(c1)              &
                        + body_fz(c2) / grid % vol(c2))             &
                        * grid % dz(s)                           
        gravity_source = ((grid % xf(s) - grid % xc(c1)) * grav_x   &
                        + (grid % yf(s) - grid % yc(c1)) * grav_y   &
                        + (grid % zf(s) - grid % zc(c1)) * grav_z)  &
                        * density(c1)                               &                                         
                        -((grid % xf(s) - grid % xc(c2)) * grav_x   &
                        + (grid % yf(s) - grid % yc(c2)) * grav_y   &
                        + (grid % zf(s) - grid % zc(c2)) * grav_z)  &
                        * density(c2)         

        a12 = 0.5 * a % fc(s)                    &
             * ( grid % vol(c1) / a % sav(c1)    &
               + grid % vol(c2) / a % sav(c2) )

        gravity_source =  a12 * (gravity_source - dotprod)           

        m_flux % n(s) = m_flux % n(s) + gravity_source

        b(c1) = b(c1) - gravity_source           
        b(c2) = b(c2) + gravity_source

      end if

    end do

  end if

  end subroutine
