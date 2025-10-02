!==============================================================================!
  subroutine Correct_Beta(Vof)
!------------------------------------------------------------------------------!
!   Step 2 of CICSAM: Correct beta for computation of volume fraction          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: fun
  type(Face_Type),  pointer :: v_flux
  real, contiguous, pointer :: beta_f(:)
  real, contiguous, pointer :: c_d(:)
  integer                   :: s, c1, c2, donor, accept
  real                      :: e_plus, e_minus, cf, delta_alfa, bcorr
!==============================================================================!

  ! Take aliases
  Flow   => Vof % pnt_flow
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  fun    => Vof % fun
  beta_f => Vof % beta_f
  c_d    => Vof % c_d

  ! Interior faces
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

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
      delta_alfa = 0.5 * (fun % n(accept) + fun % o(accept)      &
                       - (fun % n(donor)  + fun % o(donor)))

      cf = c_d(donor)

      if(fun % n(donor) < 0.0) then
        e_minus = max(-fun % n(donor), 0.0)
        ! Donor value < 0.0 Ex: sd = -0.1 -> e_minus = +0.1
        if(e_minus > FEMTO .and. cf > FEMTO) then
          if(delta_alfa > e_minus) then
            bcorr = e_minus * (2.0 + cf - 2.0 * cf * beta_f(s))     &
                            / (2.0 * cf * (delta_alfa - e_minus))
            bcorr = min(bcorr, beta_f(s))
          end if
        end if
      end if

      if(fun % n(donor) > 1.0) then
        e_plus = max(fun % n(donor) - 1.0, 0.0)
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

  end do

  end subroutine
