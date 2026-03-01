!==============================================================================!
  subroutine Vis_T_K_Omega_Sst(Turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for k-omega SST model.                    !
!                                                                              !
!   In the domain:                                                             !
!     mu_t = rho * a1 * k / max(a1*omega, S*f2)                                !
!                                                                              !
!   On the boundary (wall viscosity):                                          !
!     keep the existing wall-viscosity model based on y+ (as in k-eps).        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: kin, omega
  integer                   :: c
  real                      :: kin_vis, denom, f2, arg2
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  kin   => Turb % kin
  omega => Turb % omega

  !-----------------------------!
  !   In the domain: mu_t SST   !
  !-----------------------------!
  do c = Cells_In_Domain()

    kin_vis = Flow % viscosity(c) / Flow % density(c)

    arg2 = max( 2.0 * sqrt(max(kin % n(c), 0.0))    /         &
                max(Turb % beta_star * omega % n(c)           &
                * Grid % wall_dist(c), TINY),                 &
                500.0 * kin_vis / max(Grid % wall_dist(c)**2  &
                * omega % n(c), TINY) )

    f2 = tanh(arg2*arg2)

    denom = max(Turb % a1 * omega % n(c), Flow % shear(c) * f2)

    Turb % vis_t(c) = Flow % density(c) * Turb % a1           &
                      * kin % n(c) / max(denom, TINY)

  end do

  !-------------------!
  !   Wall function   !
  !-------------------+
  call Turb % Wall_Function()

  call Grid % Exchange_Cells_Real(Turb % vis_w)
  if(Flow % n_scalars > 0) then
    call Grid % Exchange_Cells_Real(Turb % diff_w)
  end if
  if(Flow % heat_transfer) then
    call Grid % Exchange_Cells_Real(Turb % con_w)
  end if

  call Grid % Exchange_Cells_Real(Turb % vis_t)

  end subroutine
