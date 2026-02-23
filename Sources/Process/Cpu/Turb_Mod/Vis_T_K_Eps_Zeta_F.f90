!==============================================================================!
  subroutine Vis_T_K_Eps_Zeta_F(Turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent (viscosity/density) for RANS models.                !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Turb_Type), target :: Turb
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  integer                   :: c, c1, c2, s, reg
  real                      :: u_tan, u_tau
  real                      :: beta, pr, sc
  real                      :: u_plus, ebf, kin_vis
  real                      :: z_o
!------------------------------------------------------------------------------!
!   Dimensions:                                                                !
!                                                                              !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   Wall visc.    vis_w    [kg/(m*s)]  | kinematic viscosity       [m^2/s]     !
!   Thermal cap.  capacity[m^2/(s^2*K)]| Therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum    (u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)

  call Turb % Time_And_Length_Scale(Grid)

  ! Pure k-eps-zeta-f
  do reg = Boundary_And_Inside_Regions()
    do c = Cells_In_Region(reg)
      Turb % vis_t(c) = Turb % c_mu_d * Flow % density(c) * zeta % n(c)  &
                      * kin % n(c) * Turb % t_scale(c)
    end do
  end do

  ! Hybrid between k-eps-zeta-f and dynamic SGS model
  if(Turb % model .eq. HYBRID_LES_RANS) then

    do reg = Boundary_And_Inside_Regions()
      do c = Cells_In_Region(reg)
        Turb % vis_t_eff(c) = max(Turb % vis_t(c),  &
                                  Turb % vis_t_sgs(c))
      end do
    end do
    call Grid % Exchange_Cells_Real(Turb % vis_t_eff)

  end if

  !-------------------!                                                          
  !   Wall function   !                                                          
  !-------------------+                                                          
  call Turb % Wall_Function()                                                       
                                                                                    
  call Grid % Exchange_Cells_Real(Turb % vis_w)                                     
  if(Flow % n_scalars > 0) call Grid % Exchange_Cells_Real(Turb % diff_w)           
  if(Flow % heat_transfer) call Grid % Exchange_Cells_Real(Turb % con_w)            
                                                                                    
  call Grid % Exchange_Cells_Real(Turb % vis_t) 

  end subroutine
