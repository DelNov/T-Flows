!==============================================================================!
  subroutine Multiphase_Mod_Vof_Predict_Beta(mult, beta_f, beta_c, c_d)
!------------------------------------------------------------------------------!
!   Step 1 of CICSAM: Compute beta for computation of volume fraction          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: beta_f(mult % pnt_grid % n_faces)
  real                          :: beta_c(mult % pnt_grid % n_faces)
  real                          :: c_d(-mult % pnt_grid % n_bnd_cells  &
                                       :mult % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: vof
  type(Face_Type),  pointer :: v_flux
  integer                   :: s
  integer                   :: c1, c2, donor, accept
  real                      :: fs, dotprod
  real                      :: alfa_u, alfa_d, alfa_a, alfa_d_til, alfa_cbc
  real                      :: alfa_uq, gamma_f, alfa_f_til, signo
  real                      :: alfa_superbee, alfa_stoic
  real                      :: cod, prodmag, ang, epsloc
!==============================================================================!

  ! Take aliases
  flow   => mult % pnt_flow
  grid   => flow % pnt_grid
  vof    => mult % vof
  v_flux => flow % v_flux

  epsloc = epsilon(epsloc)

  if(vof % adv_scheme .eq. CICSAM) then

    !--------------------!
    !   Compute beta_f   !
    !--------------------!

    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      beta_f(s) = 0.0
      if(abs(v_flux % n(s)) > epsloc) then

        if(v_flux % n(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        alfa_d = vof % n(donor)
        alfa_a = vof % n(accept)

        dotprod = signo * (vof % x(donor) * grid % dx(s)        &
                          + vof % y(donor) * grid % dy(s)        &
                          + vof % z(donor) * grid % dz(s))

        alfa_u = min(max(alfa_a - 2.0 * dotprod, 0.0), 1.0) !old way

!        call Multiphase_Mod_Vof_Find_Upstream_Phi(vof,      &
!                                                  vof % x,  &
!                                                  vof % y,  &
!                                                  vof % z,  &
!                                                  s, donor, accept, alfa_u)
        ! Face is inside the domain
        if(abs(alfa_u - alfa_a) > epsloc) then

          alfa_d_til = (alfa_d - alfa_u) / (alfa_a - alfa_u)

          cod = min(1.0, c_d(donor))

          ! Compute alfa_cbc
          if(alfa_d_til >= 0.0 .and. alfa_d_til <= 1.0) then
            alfa_cbc = min(1.0, alfa_d_til / max(cod, epsloc))
          else
            alfa_cbc = alfa_d_til
          end if

          ! Compute alfa_uq
          if(alfa_d_til >= 0.0 .and. alfa_d_til <= 1.0) then
            alfa_uq = min(cod                           &
                        * alfa_d_til + (1.0 - cod)      &
                        * (6.0 * alfa_d_til + 3.0)      &
                         / 8.0, alfa_cbc)
          else
            alfa_uq = alfa_d_til
          end if

          ! Compute angle:
          prodmag = sqrt(vof % x(donor) ** 2              &
                       + vof % y(donor) ** 2              &
                       + vof % z(donor) ** 2)             &
                  * grid % d(s)

          if(prodmag > epsloc) then
            ang = dotprod / prodmag
          else
            ang = 1.0 / epsloc
          end if

          gamma_f = min(ang ** 2.0, 1.0)

          alfa_f_til = gamma_f * alfa_cbc + (1.0 - gamma_f) * alfa_uq

          if(abs(1.0 - alfa_d_til) > epsloc) then

            beta_f(s) = min(max((alfa_f_til - alfa_d_til)                   &
                              / (1.0 - alfa_d_til) + beta_c(s), 0.0), 1.0)

          end if
        end if
      end if
    end do

  else if(vof % adv_scheme .eq. STACS) then

    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      beta_f(s) = 0.0
      if(abs(v_flux % n(s)) > epsloc) then

        if(v_flux % n(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        alfa_d = vof % n(donor)
        alfa_a = vof % n(accept)

        dotprod = signo * (vof % x(donor) * grid % dx(s)        &
                          + vof % y(donor) * grid % dy(s)        &
                          + vof % z(donor) * grid % dz(s))

        ! Face is inside the domain
        alfa_u = min(max(alfa_a - 2.0 * dotprod, 0.0), 1.0)

        if(abs(alfa_u - alfa_a) > epsloc) then

          alfa_d_til = (alfa_d - alfa_u) / (alfa_a - alfa_u)

          ! Compute alfa_SUBERBEE
          if(alfa_d_til > 0.0 .and. alfa_d_til < 1.0) then
            alfa_superbee = 1.0
          else
            alfa_superbee = alfa_d_til
          end if

          if(alfa_d_til > 0.0 .and. alfa_d_til <= 0.5) then
            alfa_stoic = 0.5 + 0.5 * alfa_d_til
          else if(alfa_d_til > 0.5 .and. alfa_d_til <= 5.0 / 6.0) then
            alfa_stoic = 3.0 / 8.0 + 3.0 / 4.0 * alfa_d_til
          else if(alfa_d_til > 5.0 / 6.0 .and. alfa_d_til <= 1.0) then
            alfa_stoic = 1.0
          else
            alfa_stoic = alfa_d_til
          end if

          ! Compute angle:
          prodmag = sqrt(vof % x(donor) ** 2              &
                       + vof % y(donor) ** 2              &
                       + vof % z(donor) ** 2)             &
                  * grid % d(s)

          if(prodmag > epsloc) then
            ang = dotprod / prodmag
          else
            ang = 1.0 / epsloc
          end if

          gamma_f = min(ang ** 4.0, 1.0)

          alfa_f_til = gamma_f * alfa_superbee + (1.0 - gamma_f) * alfa_stoic

          if(abs(1.0 - alfa_d_til) > epsloc) then

            beta_f(s) = min(max((alfa_f_til - alfa_d_til)                   &
                              / (1.0 - alfa_d_til) + beta_c(s), 0.0), 1.0)

          end if

        end if
      end if
    end do
  end if

  end subroutine
