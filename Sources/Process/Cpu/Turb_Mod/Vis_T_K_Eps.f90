!==============================================================================!
  subroutine Vis_T_K_Eps(Turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!                                                                              !
!   In the domain:                                                             !
!   For k-eps model :                                                          !
!                                                                              !
!   vis_t = c_mu * rho * k^2 * eps                                             !
!                                                                              !
!   On the boundary (wall viscosity):                                          !
!   vis_tw = y^+ * vis_t kappa / (E * ln(y^+))                                 !
!                                                                              !
!   For k-eps-v2f model :                                                      !
!                                                                              !
!   vis_t = CmuD * rho * Tsc  * vv                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: kin, eps
  integer                   :: c1, c2, s, c, reg
  real                      :: pr, beta, ebf, sc
  real                      :: u_tan, u_tau
  real                      :: kin_vis, u_plus, y_star, re_t, f_mu
  real                      :: z_o
!------------------------------------------------------------------------------!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear    [1/s]        !
!   dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t    [kg/(m*s)]   !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity[kg/(m*s)]   !
!   density       density  [kg/m^3]    | Turb. kin en.   kin % n  [m^2/s^2]    !
!   cell volume   vol      [m^3]       | length          lf       [m]          !
!   left hand s.  A        [kg/s]      | right hand s.   b        [kg*m^2/s^3] !
!   wall visc.    vis_w    [kg/(m*s)]  | kinematic viscosity      [m^2/s]      !
!   thermal cap.  capacity[m^2/(s^2*K)]| therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  call Turb % Alias_K_Eps   (kin, eps)

  do c = Cells_In_Domain()

    ! Kinematic viscosities
    kin_vis = Flow % viscosity(c) / Flow % density(c)

    re_t =  Flow % density(c) * kin % n(c)**2  &
         / (Flow % viscosity(c) * eps % n(c))

    y_star = (kin_vis * eps % n(c))**0.25 * Grid % wall_dist(c)/kin_vis

    f_mu = (1.0 -     exp(-y_star/14.0))**2   &
         * (1.0 + 5.0*exp(-(re_t/200.0) * (re_t/200.0) ) /re_t**0.75)

    f_mu = min(1.0,f_mu)

    Turb % vis_t(c) = f_mu * Turb % c_mu * Flow % density(c) * kin % n(c)**2  &
                      / (eps % n(c) + TINY)
  end do

  !-------------------!
  !   Wall function   !
  !-------------------+
  call Turb % Wall_Function()

  call Grid % Exchange_Cells_Real(Turb % vis_w)
  if(Flow % heat_transfer) then
    call Grid % Exchange_Cells_Real(Turb % con_w)
  end if
  if(Flow % n_scalars > 0) then
    call Grid % Exchange_Cells_Real(Turb % diff_w)
  end if

  call Grid % Exchange_Cells_Real(Turb % vis_t)

  end subroutine
