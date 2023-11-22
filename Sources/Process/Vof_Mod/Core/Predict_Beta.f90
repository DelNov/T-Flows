!==============================================================================!
  subroutine Predict_Beta(Vof)
!------------------------------------------------------------------------------!
!   Step 1 of CICSAM: Compute beta for computation of volume fraction          !
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
  real, contiguous, pointer :: beta_c(:)
  real, contiguous, pointer :: c_d(:)
  integer                   :: s
  integer                   :: c1, c2, donor, accept
  real                      :: dotprod
  real                      :: alfa_u, alfa_d, alfa_a, alfa_d_til, alfa_cbc
  real                      :: alfa_uq, gamma_f, alfa_f_til, signo
  real                      :: alfa_superbee, alfa_stoic
  real                      :: cod, prodmag, ang
!==============================================================================!

  ! Take aliases
  Flow   => Vof % pnt_flow
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  fun    => Vof % fun
  beta_f => Vof % beta_f
  beta_c => Vof % beta_c
  c_d    => Vof % c_d

  if(fun % adv_scheme .eq. CICSAM) then

    !--------------------!
    !   Compute beta_f   !
    !--------------------!

    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      beta_f(s) = 0.0
      if(abs(v_flux % n(s)) > FEMTO) then

        if(v_flux % n(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        alfa_d = fun % n(donor)
        alfa_a = fun % n(accept)

        dotprod = signo * (  fun % x(donor) * Grid % dx(s)   &
                           + fun % y(donor) * Grid % dy(s)   &
                           + fun % z(donor) * Grid % dz(s))

        alfa_u = min(max(alfa_a - 2.0 * dotprod, 0.0), 1.0)  ! old way

!       call Multiphase_Mod_Vof_Find_Upstream_Phi(fun,      &
!                                                 fun % x,  &
!                                                 fun % y,  &
!                                                 fun % z,  &
!                                                 s, donor, accept, alfa_u)
        ! Face is inside the domain
        if(abs(alfa_u - alfa_a) > FEMTO) then

          alfa_d_til = (alfa_d - alfa_u) / (alfa_a - alfa_u)

          cod = min(1.0, c_d(donor))

          ! Compute alfa_cbc
          if(alfa_d_til >= 0.0 .and. alfa_d_til <= 1.0) then
            alfa_cbc = min(1.0, alfa_d_til / max(cod, FEMTO))
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
          prodmag = sqrt(  fun % x(donor) ** 2   &
                         + fun % y(donor) ** 2   &
                         + fun % z(donor) ** 2)  &
                  * Grid % d(s)

          if(prodmag > FEMTO) then
            ang = dotprod / prodmag
          else
            ang = 1.0 / FEMTO
          end if

          gamma_f = min(ang ** 2, 1.0)

          alfa_f_til = gamma_f * alfa_cbc + (1.0 - gamma_f) * alfa_uq

          if(abs(1.0 - alfa_d_til) > FEMTO) then

            beta_f(s) = min(max((alfa_f_til - alfa_d_til)                   &
                              / (1.0 - alfa_d_til) + beta_c(s), 0.0), 1.0)

          end if
        end if
      end if    ! abs(v_flux % n(s)) > FEMTO
    end do      ! faces

  else if(fun % adv_scheme .eq. STACS) then

    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      beta_f(s) = 0.0
      if(abs(v_flux % n(s)) > FEMTO) then

        if(v_flux % n(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        alfa_d = fun % n(donor)
        alfa_a = fun % n(accept)

        dotprod = signo * (  fun % x(donor) * Grid % dx(s)   &
                           + fun % y(donor) * Grid % dy(s)   &
                           + fun % z(donor) * Grid % dz(s))

        ! Face is inside the domain
        alfa_u = min(max(alfa_a - 2.0 * dotprod, 0.0), 1.0)

        if(abs(alfa_u - alfa_a) > FEMTO) then

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
          prodmag = sqrt(  fun % x(donor) ** 2   &
                         + fun % y(donor) ** 2   &
                         + fun % z(donor) ** 2)  &
                  * Grid % d(s)

          if(prodmag > FEMTO) then
            ang = dotprod / prodmag
          else
            ang = 1.0 / FEMTO
          end if

          gamma_f = min(ang ** 4, 1.0)

          alfa_f_til = gamma_f * alfa_superbee + (1.0 - gamma_f) * alfa_stoic

          if(abs(1.0 - alfa_d_til) > FEMTO) then

            beta_f(s) = min(max((alfa_f_til - alfa_d_til)                   &
                              / (1.0 - alfa_d_til) + beta_c(s), 0.0), 1.0)

          end if

        end if

      end if  ! abs(v_flux % n(s)) > FEMTO
    end do    ! faces

  end if      ! STACS

  end subroutine
