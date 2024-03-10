!==============================================================================!
  subroutine Vis_T_Rsm(Turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RSM models ('EBM' and 'HJ').          !
!   If hybrid option is used turbulent diffusivity is modeled by vis_t.        !
!   Otherwise, vis_t is used as false diffusion in order to increase           !
!   stability of computation.                                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: kin, eps
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  integer                   :: c
  real                      :: cmu_mod
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   wall visc.    vis_w    [kg/(m*s)]  | kinematic viscosity       [m^2/s]     !
!   thermal cap.  capacity[m^2/(s^2*K)]| therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  call Turb % Alias_K_Eps   (kin, eps)
  call Turb % Alias_Stresses(uu, vv, ww, uv, uw, vw)

  call Calculate_Shear_And_Vorticity(Flow)

  do c = Cells_In_Domain()
    kin % n(c) = 0.5*max(uu % n(c) + vv % n(c) + ww % n(c), TINY)

    cmu_mod = max(-(  uu % n(c) * u % x(c)               &
                    + vv % n(c) * v % y(c)               &
                    + ww % n(c) * w % z(c)               &
                    + uv % n(c) * (v % x(c) + u % y(c))  &
                    + uw % n(c) * (u % z(c) + w % x(c))  &
                    + vw % n(c) * (v % z(c) + w % y(c))) &
      / (kin % n(c) * Turb % t_scale(c) * Flow % shear(c)**2 + TINY), 0.0)

    cmu_mod = min(0.12, cmu_mod)
    Turb % vis_t(c) = cmu_mod * Flow % density(c)  &
                    * kin % n(c) * Turb % t_scale(c)
    Turb % vis_t(c) = max(Turb % vis_t(c), TINY)
  end do

  call Grid % Exchange_Cells_Real(Turb % vis_t)

  end subroutine
