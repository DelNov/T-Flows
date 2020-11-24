!==============================================================================!
  subroutine Multiphase_Mod_Vof_Correct_Beta(mult, beta_f, c_d)
!------------------------------------------------------------------------------!
!   Step 2 of CICSAM: Correct beta for computation of volume fraction          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: beta_f(mult % pnt_grid % n_faces)
  real                          :: c_d(-mult % pnt_grid % n_bnd_cells  &
                                       :mult % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: vof
  type(Face_Type),  pointer :: v_flux
  integer                   :: s, c1, c2, donor, accept
  real                      :: e_plus, e_minus, cf, delta_alfa, bcorr
!==============================================================================!

  ! Take aliases
  flow   => mult % pnt_flow
  grid   => flow % pnt_grid
  vof    => mult % vof
  v_flux => flow % v_flux

  ! Interior faces
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 > 0) then

      if(abs(v_flux % n(s)) > FEMTO) then

        if(v_flux % n(s) > 0.0) then
          donor = c1
          accept = c2
        else
          donor = c2
          accept = c1
        end if

        !--------------------!
        !   Correct beta_f   !
        !--------------------!
        bcorr = 0.0
        delta_alfa = 0.5 * (vof % n(accept) + vof % o(accept)      &
                         - (vof % n(donor)  + vof % o(donor)))

        cf = c_d(donor)

        if(vof % n(donor) < 0.0) then
          e_minus = max(-vof % n(donor), 0.0)
          ! Donor value < 0.0 Ex: sd = -0.1 -> e_minus = +0.1
          if(e_minus > FEMTO .and. cf > FEMTO) then
            if(delta_alfa > e_minus) then
              bcorr = e_minus * (2.0 + cf - 2.0 * cf * beta_f(s))     &
                              / (2.0 * cf * (delta_alfa - e_minus))
              bcorr = min(bcorr, beta_f(s))
            end if
          end if
        end if

        if(vof % n(donor) > 1.0) then
          e_plus = max(vof % n(donor) - 1.0, 0.0)
          ! Donor value > 1.0 Ex: sd = 1.1 -> e_plus = +0.1
          if(e_plus > FEMTO .and. cf > FEMTO) then
            if(delta_alfa < - e_plus) then
              bcorr = e_plus * (2.0 + cf - 2.0 * cf * beta_f(s))     &
                             / (2.0 * cf * (-delta_alfa - e_plus))
              bcorr = min(bcorr, beta_f(s))
            end if
          end if

        end if

        beta_f(s) = beta_f(s) - bcorr
        beta_f(s) = max(beta_f(s), 0.0)

      end if

    end if  ! c2 > 0

  end do

  end subroutine
