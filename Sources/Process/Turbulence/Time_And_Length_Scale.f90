!==============================================================================!
  subroutine Time_And_Length_Scale(grid)
!------------------------------------------------------------------------------!
!   Calculates time scale and leght scale in manner to avoid singularity       !
!   in eps equation.                                                           !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Control_Mod
  use Work_Mod, only: t1    => r_cell_01,  &  ! [s]
                      t2    => r_cell_02,  &  ! [s]
                      t3    => r_cell_03,  &  ! [s]
                      l1    => r_cell_04,  &  ! [m]
                      l2    => r_cell_05,  &  ! [m]
                      l3    => r_cell_06,  &  ! [m]
                      eps_l => r_cell_07      ! [m]
!------------------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!----------------------------------[Locals]------------------------------------!
  real    :: kin_vis   ! kinematic viscosity [m^2/s]
  integer :: c
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!

  kin_vis = viscosity/density
  eps_l(1:grid % n_cells) = eps % n(1:grid % n_cells) + TINY ! limited eps % n

  do c = 1, grid % n_cells
    t1(c) = kin % n(c)/eps_l(c)
    t2(c) = c_t*sqrt(kin_vis/eps_l(c))

    l1(c) = kin % n(c)**1.5/eps_l(c)
    l2(c) = c_nu * (kin_vis**3 / eps_l(c))**0.25
  end do

  if(turbulence_model .eq. K_EPS_ZETA_F) then

    if(rough_walls) then
      do c = 1, grid % n_cells
        t_scale(c) =     max(t1(c),t2(c))
        l_scale(c) = c_l*max(l1(c),l2(c))
      end do
    else
      do c = 1, grid % n_cells
        t3(c) = 0.6/(sqrt(3.0)*c_mu_d * zeta % n(c) * shear(c) + TINY)
        l3(c) = sqrt(kin % n(c)/3.0)/(c_mu_d * zeta % n(c) * shear(c) + TINY)
      end do
      do c = 1, grid % n_cells
        t_scale(c) =     max(min(t1(c),t3(c)),t2(c))
        l_scale(c) = c_l*max(min(l1(c),l3(c)),l2(c))
      end do
    end if

  else if(turbulence_model .eq. REYNOLDS_STRESS) then
    do c = 1, grid % n_cells
      kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)
      t_scale(c) =     max(t1(c),t2(c))
      l_scale(c) = c_l*max(l1(c),l2(c))
    end do
  end if

  end subroutine
