!==============================================================================!
  subroutine Multiphase_Mod_Vof_Skewness_Correction(mult, grid, beta_c)
!------------------------------------------------------------------------------!
!   Implicit skewness correction for beta_f, based on the work of Denner 2014  !
!   https://doi.org/10.1016/j.jcp.2014.09.002
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Grid_Type)               :: grid
  real                          :: beta_c(1:grid % n_faces)  !Skew. correction
  real                          :: c_d(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),     pointer :: flow
  type(Var_Type),       pointer :: vof
  type(Face_Type),      pointer :: m_flux
  integer                       :: s
  integer                       :: c1, c2, donor, accept
  real                          :: fs, dot_prod
  real                          :: alfa_u, alfa_d, alfa_a, alfa_d_til, alfa_cbc
  real                          :: alfa_uq, gamma_f, alfa_f_til, signo
  real                          :: alfa_superbee, alfa_stoic
  real                          :: cod, prodmag, ang, epsloc
  real                          :: gf_x, gf_y, gf_z  ! avg gradient
  real                          :: gradf_x, gradf_y, gradf_z  ! gradient at f
  real                          :: corr_x, corr_y, corr_z
  real                          :: sd, nplane(3), uline(3)
  real                          :: r_f(3), f_prime(3)
  real                          :: dsc1, dsc2, dotprod1, dotprod2
  !==============================================================================!

  ! Take aliases
  flow   => mult % pnt_flow
  vof    => mult % vof
  m_flux => flow % m_flux

  epsloc = epsilon(epsloc)

  !--------------------!
  !   Compute beta_c   !
  !--------------------!

  beta_c = 0.0

  if (mult % skew_corr) then
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      if (abs(m_flux % n(s)) > epsloc) then

        if (m_flux % n(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        if (abs(vof % n(c1) - vof % n(c2)) > epsloc) then
          ! Find intersection point on plane (f')
          ! https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
          nplane = (/grid % sx(s), grid % sy(s), grid % sz(s)/)
          uline = (/grid % dx(s), grid % dy(s), grid % dz(s)/)

          sd = dot_product(nplane, (/grid % xf(s) - grid % xc(c1),     &
                                     grid % yf(s) - grid % yc(c1),     &
                                     grid % zf(s) - grid % zc(c1)/) )  &
                                     /dot_product(nplane,uline)

          f_prime = (/grid % xc(c1) + sd * uline(1),  &
                      grid % yc(c1) + sd * uline(2),  &
                      grid % zc(c1) + sd * uline(3)/)

          dsc1 = Math_Mod_Distance(grid % xc(c1), grid % yc(c1), grid % zc(c1),   &
                                   f_prime(1), f_prime(2), f_prime(3))

          dsc2 = Math_Mod_Distance(grid % xc(c2), grid % yc(c2), grid % zc(c2),   &
                                   f_prime(1), f_prime(2), f_prime(3))

          ! Interpolation factor
          !fs = dsc2 / (dsc1 + dsc2)
          fs = grid % f(s)
          r_f = (/grid % xf(s) - f_prime(1),  &
                  grid % yf(s) - f_prime(2),  &
                  grid % zf(s) - f_prime(3)/)

          if (norm2(r_f) > epsloc) then
            ! Average Gradient
            gf_x = fs * vof % x(c1) + (1.0 - fs) * vof % x(c2)
            gf_y = fs * vof % y(c1) + (1.0 - fs) * vof % y(c2)
            gf_z = fs * vof % z(c1) + (1.0 - fs) * vof % z(c2)

            dotprod1 = dot_product((/grid % sx(s), grid % sy(s), grid % sz(s)/),  &
                                   (/grid % dx(s), grid % dy(s), grid % dz(s)/))

            dotprod2 = dot_product((/gf_x, gf_y, gf_z/),  &
                                   (/grid % dx(s), grid % dy(s), grid % dz(s)/))

            gradf_x = gf_x + grid % sx(s) / dotprod1 * (vof % n(c2) - vof % n(c1)  &
                                                       - dotprod2)

            gradf_y = gf_y + grid % sy(s) / dotprod1 * (vof % n(c2) - vof % n(c1)  &
                                                       - dotprod2)

            gradf_z = gf_z + grid % sz(s) / dotprod1 * (vof % n(c2) - vof % n(c1)  &
                                                       - dotprod2)

            beta_c(s) = dot_product((/gradf_x, gradf_y, gradf_z/),r_f) &
                      / (vof % n(accept) - vof % n(donor))
          end if
        end if
      end if
    end do
  end if

  end subroutine
