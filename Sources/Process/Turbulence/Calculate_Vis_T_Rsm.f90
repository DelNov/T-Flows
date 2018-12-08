!==============================================================================!
  subroutine Calculate_Vis_T_Rsm(flow)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RSM models ('EBM' and 'HJ').          !
!   If hybrid option is used turbulent diffusivity is modeled by vis_t.        !
!   Otherwise, vis_t is used as false diffusion in order to increase           !
!   stability of computation.                                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: c
  real                     :: cmu_mod
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   wall visc.    vis_wall [kg/(m*s)]  | kinematic viscosity       [m^2/s]     !
!   thermal cap.  capacity[m^2/(s^2*K)]| therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w

  call Calculate_shear_And_Vorticity(flow)

  if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    do c = 1, grid % n_cells
      kin % n(c) = 0.5*max(uu % n(c) + vv % n(c) + ww % n(c), TINY)

      cmu_mod = max(-(  uu % n(c) * u % x(c)               &
                      + vv % n(c) * v % y(c)               &
                      + ww % n(c) * w % z(c)               &
                      + uv % n(c) * (v % x(c) + u % y(c))  &
                      + uw % n(c) * (u % z(c) + w % x(c))  &
                      + vw % n(c) * (v % z(c) + w % y(c))) &
        / max(kin % n(c)**2 / max(eps_tot(c), TINY) * shear(c)**2, TINY), 0.0)

      cmu_mod = min(0.12, cmu_mod)
      vis_t(c) = cmu_mod * density * kin % n(c)**2 / max(eps_tot(c), TINY)
      vis_t(c) = max(vis_t(c), TINY)
    end do
  else if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
    do c = 1, grid % n_cells
      kin % n(c) = 0.5*max(uu % n(c) + vv % n(c) + ww % n(c), TINY)

      ! Pk/ (k^2/eps * S^2)
      cmu_mod = max(-(  uu % n(c) * u % x(c)                         &
                      + vv % n(c) * v % y(c)                         &
                      + ww % n(c) * w % z(c)                         &
                      + uv % n(c) * (v % x(c) + u % y(c))            &
                      + uw % n(c) * (u % z(c) + w % x(c))            &
                      + vw % n(c) * (v % z(c) + w % y(c)))           &
                      / max(  kin % n(c)**2 / max(eps % n(c), TINY)  &
                              * shear(c)**2, TINY), 0.0)

      cmu_mod = min(0.12,cmu_mod)
      vis_t(c) = cmu_mod * density * kin % n(c)**2 / max(eps % n(c), TINY)
      vis_t(c) = max(vis_t(c), TINY)
    end do
  end if

  call Comm_Mod_Exchange_Real(grid, vis_t)

  end subroutine
