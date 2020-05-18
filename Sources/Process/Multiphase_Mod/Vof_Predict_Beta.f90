!==============================================================================!
  subroutine Multiphase_Mod_Vof_Predict_Beta(phi,                  &
                                             phi_flux,             &
                                             phi_x, phi_y, phi_z,  &
                                             di, dj, dk,           &
                                             beta_f,               &
                                             c_d)
!------------------------------------------------------------------------------!
!   Step 1 of CICSAM: Compute beta for computation of volume fraction          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
  real           :: phi_flux(phi % pnt_grid % n_faces)
  real           :: phi_x(-phi % pnt_grid % n_bnd_cells:  &
                           phi % pnt_grid % n_cells),     &
                    phi_y(-phi % pnt_grid % n_bnd_cells:  &
                           phi % pnt_grid % n_cells),     &
                    phi_z(-phi % pnt_grid % n_bnd_cells:  &
                           phi % pnt_grid % n_cells)
  real           :: di(phi % pnt_grid % n_faces),         &
                    dj(phi % pnt_grid % n_faces),         &
                    dk(phi % pnt_grid % n_faces)
  real           :: beta_f(phi % pnt_grid % n_faces)
  real           :: c_d(-phi % pnt_grid % n_bnd_cells:  &
                         phi % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s
  integer                  :: c1, c2, donor, accept
  real                     :: fs, dot_prod
  real                     :: alfa_u, alfa_d, alfa_a, alfa_d_til, alfa_cbc
  real                     :: alfa_uq, gamma_f, alfa_f_til, signo
  real                     :: alfa_superbee, alfa_stoic
  real                     :: cod, prodmag, ang, epsloc
!==============================================================================!

  ! Take aliases
  grid => phi % pnt_grid

  epsloc = epsilon(epsloc)

  if(phi % adv_scheme .eq. CICSAM) then

    !--------------------!
    !   Compute beta_f   !
    !--------------------!

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      beta_f(s) = 0.0
      if (abs(phi_flux(s)) > epsloc) then 

        if (phi_flux(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        if (c2 > 0) then
          alfa_d = phi % n(donor)
          alfa_a = phi % n(accept)

          dot_prod = signo * (phi_x(donor) * di(s)        &
                            + phi_y(donor) * dj(s)        &
                            + phi_z(donor) * dk(s))

          ! Face is inside the domain
          alfa_u = min(max(alfa_a - 2.0 * dot_prod, 0.0), 1.0) !old way

!          call Multiphase_Mod_Vof_Find_Upstream_phi(phi,    &
!                                                    phi_x,  &
!                                                    phi_y,  &
!                                                    phi_z,  &    
!                                                    s, donor, accept, alfa_u)
          if (abs(alfa_u - alfa_a) > TINY) then

            alfa_d_til = (alfa_d - alfa_u) / (alfa_a - alfa_u)

            cod = min(1.0, c_d(donor))

            ! Compute alfa_cbc 
            if (alfa_d_til >= 0.0 .and. alfa_d_til <= 1.0) then
              alfa_cbc = min(1.0, alfa_d_til / cod)
            else
              alfa_cbc = alfa_d_til
            end if

            ! Compute alfa_uq
            if (alfa_d_til >= 0.0 .and. alfa_d_til <= 1.0) then
              alfa_uq = min(cod                           & 
                          * alfa_d_til + (1.0 - cod)      &
                          * (6.0 * alfa_d_til + 3.0)      &
                           / 8.0, alfa_cbc)
            else
              alfa_uq = alfa_d_til
            end if

            ! Compute angle:
            prodmag = sqrt(phi_x(donor) ** 2              &
                         + phi_y(donor) ** 2              &
                         + phi_z(donor) ** 2)             &
                    * sqrt(di(s) ** 2 + dj(s) ** 2 + dk(s) ** 2)

            if (prodmag > TINY) then
              ang = dot_prod / prodmag
            else
              ang = 1.0 / TINY
            end if

            gamma_f = min(ang ** 2.0, 1.0)

            alfa_f_til = gamma_f * alfa_cbc + (1.0 - gamma_f) * alfa_uq

            if (abs(1.0 - alfa_d_til) > TINY) then

              beta_f(s) = min(max((alfa_f_til - alfa_d_til)                   &
                                / (1.0 - alfa_d_til), 0.0), 1.0)

            end if
          end if
        end if
      end if
    end do

  else if (phi % adv_scheme .eq. STACS) then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      beta_f(s) = 0.0
      if (abs(phi_flux(s)) > epsloc) then 

        if (phi_flux(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        if (c2 > 0) then

          alfa_d = phi % n(donor)
          alfa_a = phi % n(accept)

          dot_prod = signo * (phi_x(donor) * di(s)        &
                            + phi_y(donor) * dj(s)        &
                            + phi_z(donor) * dk(s))

          ! Face is inside the domain
          alfa_u = min(max(alfa_a - 2.0 * dot_prod, 0.0), 1.0)

          if (abs(alfa_u - alfa_a) > TINY) then

            alfa_d_til = (alfa_d - alfa_u) / (alfa_a - alfa_u)

            ! Compute alfa_SUBERBEE
            if (alfa_d_til > 0.0 .and. alfa_d_til < 1.0) then
              alfa_superbee = 1.0
            else 
              alfa_superbee = alfa_d_til
            end if

            if (alfa_d_til > 0.0 .and. alfa_d_til <= 0.5) then
              alfa_stoic = 0.5 + 0.5 * alfa_d_til
            else if (alfa_d_til > 0.5 .and. alfa_d_til <= 5.0 / 6.0) then
              alfa_stoic = 3.0 / 8.0 + 3.0 / 4.0 * alfa_d_til
            else if (alfa_d_til > 5.0 / 6.0 .and. alfa_d_til <= 1.0) then
              alfa_stoic = 1.0
            else
              alfa_stoic = alfa_d_til
            end if

            ! Compute angle:
            prodmag = sqrt(phi_x(donor) ** 2              &
                         + phi_y(donor) ** 2              &
                         + phi_z(donor) ** 2)             &     
                    * sqrt(di(s) ** 2 + dj(s) ** 2 + dk(s) ** 2)


            if (prodmag > TINY) then
              ang = dot_prod / prodmag
            else
              ang = 1.0 / tiny
            end if

            gamma_f = min(ang ** 4.0, 1.0)

            alfa_f_til = gamma_f * alfa_superbee + (1.0 - gamma_f) * alfa_stoic

            if (abs(1.0 - alfa_d_til) > TINY) then

              beta_f(s) = min(max((alfa_f_til - alfa_d_til)                   &
                                / (1.0 - alfa_d_til), 0.0), 1.0)

            end if

          end if
        end if
      end if
    end do
  end if

  end subroutine
