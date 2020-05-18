!==============================================================================!
  subroutine Multiphase_Mod_Vof_Correct_Beta(mult,      &
                                             phi,       &
                                             phi_flux,  &
                                             beta_f,    &
                                             dt)
!------------------------------------------------------------------------------!
!   Step 2 of CICSAM: Correct beta for computation of volume fraction          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type)    :: mult
  type(Var_Type)           :: phi
  real,        allocatable :: phi_flux(:)
  real,        allocatable :: beta_f(:)
  real                     :: dt
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: s, c1, c2, donor, accept
  real                      :: fs, signo, e_plus, e_minus, cf, delta_alfa, bcorr
!==============================================================================!

  ! Take aliases
  flow => mult % pnt_flow
  grid => phi  % pnt_grid

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    if (c2 > 0) then
      if (abs(phi_flux(s)) > PICO) then

        if (phi_flux(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        !--------------------!
        !   Correct beta_f   !
        !--------------------!
        bcorr = 0.0
        delta_alfa = 0.5 * (phi % n(accept) + phi % o(accept)      &
                         - (phi % n(donor) + phi % o(donor)))

        cf = min(max(- phi_flux(s) / flow % density_f(s)           &
                              * dt / grid % vol(donor),0.0),1.0)

        if (phi % n(donor) < 0.0) then
          e_minus = max(-phi % n(donor),0.0)
          ! Donor value < 0.0 Ex: sd = -0.1 -> e_minus = +0.1
          if (e_minus > TINY .and. cf > TINY) then
            if (delta_alfa > e_minus) then
              bcorr = e_minus * (2.0 + cf - 2.0 * cf * beta_f(s))     &
                              / (2.0 * cf * (delta_alfa - e_minus))
              bcorr = min(bcorr, beta_f(s))
            end if
          end if
        end if

        if (phi % n(donor) > 1.0) then
          e_plus = max(phi % n(donor) - 1.0,0.0)
          ! Donor value > 1.0 Ex: sd = 1.1 -> e_plus = +0.1
          if (e_plus > TINY .and. cf > TINY) then
            if (delta_alfa < - e_plus) then
              bcorr = e_plus * (2.0 + cf - 2.0 * cf * beta_f(s))     &
                             / (2.0 * cf * (-delta_alfa - e_plus))
              bcorr = min(bcorr, beta_f(s))
            end if
          end if

        end if

        beta_f(s) = beta_f(s) - bcorr
        beta_f(s) = max(beta_f(s),0.0)

      end if
    end if

  end do

  end subroutine
