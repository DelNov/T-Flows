!==============================================================================!
  subroutine Calculate_Vis_T_Rsm(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RSM models ('EBM' and 'HJ').          !
!   If hybrid option is used turbulent diffusivity is modeled by vis_t.        !
!   Otherwise, vis_t is used as false diffusion in order to increase           !
!   stability of computation.                                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c
  real              :: cmu_mod                                        
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   Wall visc.    vis_wall [kg/(m*s)]  |                                       !
!   Thermal cap.  capacity[m^2/(s^2*K)]| Therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!

  call Calculate_shear_And_Vorticity(grid)

  if(turbulence_model .eq. HANJALIC_JAKIRLIC) then
    do c = 1, grid % n_cells
      kin % n(c) = 0.5*max(uu % n(c) + vv % n(c) + ww % n(c), TINY)

      cmu_mod = max(-(  uu % n(c) * u % x(c)               &
                      + vv % n(c) * v % y(c)               &
                      + ww % n(c) * w % z(c)               &
                      + uv % n(c) * (v % x(c) + u % y(c))  &
                      + uw % n(c) * (u % z(c) + w % x(c))  &
                      + vw % n(c) * (v % z(c) + w % y(c))) &
        / max(kin % n(c)**2. / ( eps_tot(c) + TINY) * shear(c)**2., TINY), 0.0)

      cmu_mod = min(0.12, cmu_mod) 
      vis_t(c) = cmu_mod * density * kin % n(c)**2. / ( eps_tot(c) + TINY)
    end do 
  else if(turbulence_model .eq. REYNOLDS_STRESS) then
    do c = 1, grid % n_cells
      kin % n(c) = 0.5*max(uu % n(c) + vv % n(c) + ww % n(c), TINY)

      ! Pk/ (k^2/eps * S^2)
      cmu_mod = max(-(  uu % n(c) * u % x(c)  &
                      + vv % n(c) * v % y(c)  &
                      + ww % n(c) * w % z(c)  &
                      + uv % n(c) * (v % x(c) + u % y(c))  &
                      + uw % n(c) * (u % z(c) + w % x(c))  &
                      + vw % n(c) * (v % z(c) + w % y(c))) &
               / max(kin % n(c)**2. / eps % n(c) * shear(c)**2., TINY), 0.0)

      cmu_mod = min(0.12,cmu_mod)
      vis_t(c) = cmu_mod * density * kin % n(c)**2. / (eps % n(c) + TINY)
    end do
  end if

  call Comm_Mod_Exchange_Real(grid, vis_t)  

  end subroutine
