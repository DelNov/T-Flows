!==============================================================================!
  subroutine Update_K_Omega_Sst_Fields(Turb)
!------------------------------------------------------------------------------!
!   Computes SST blending functions F1/F2 per cell.
!
!   This routine is designed to NOT touch wall-function / heat-transfer
!   coefficients (kappa, e_log, c_theta, afm_*, ...). It only updates SST
!   blending helpers used by diffusion/viscosity/source terms.
!
!   Implemented (robust, no-gradients version):
!     F2 from arg2 (used in viscosity limiter)
!     F1 from arg1 using the first two standard max-terms
!       (the cross-diffusion-based third term can be added if grad(k)·grad(ω)
!        is available in your code; see TODO block).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: kin, omega
  integer                   :: c
  real                      :: y, nu, omg, k, arg1, arg2, t1, t2
  real, parameter           :: TINY = 1.0e-12
!==============================================================================!

  Flow  => Turb % pnt_flow
  Grid  => Flow % pnt_grid
  kin   => Turb % kin
  omega => Turb % omega

  if(.not. allocated(Turb % sst_f1)) then
    allocate(Turb % sst_f1(-Grid % n_bnd_cells:Grid % n_cells))
    allocate(Turb % sst_f2(-Grid % n_bnd_cells:Grid % n_cells))
  end if

  do c = 1, Grid % n_cells

    y   = max(Grid % wall_dist(c), TINY)
    nu  = Flow % viscosity(c) / max(Flow % density(c), TINY)
    k   = max(kin % n(c), 0.0)
    omg = max(omega % n(c), TINY)

    !-----------------------------!
    !   F2 (Menter SST: arg2)     !
    !-----------------------------!
    t1   = 2.0 * sqrt(k) / (Turb % beta_star * omg * y)
    t2   = 500.0 * nu / (y*y * omg)
    arg2 = max(t1, t2)
    Turb % sst_f2(c) = tanh(arg2*arg2)

    !-----------------------------!
    !   F1 (Menter SST: arg1)     !
    !-----------------------------!
    ! Standard definition uses:
    ! arg1 = min( max( t1, t2, t3 ), 10 * nu / (y^2 * omg) )
    ! Here we use a robust version without t3 (cross-diffusion) term:
    !   t3 = 4 * rho * sigma_w2 * k / (CD_kw * y^2)
    !
    ! TODO (optional): if you can compute CD_kw, insert t3 into max().
    arg1 = max(t1, t2)
    arg1 = min(arg1, 10.0 * nu / (y*y * omg))
    Turb % sst_f1(c) = tanh(arg1*arg1*arg1*arg1)

  end do

  end subroutine Update_K_Omega_Sst_Fields
